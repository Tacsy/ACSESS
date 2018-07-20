#!/usr/bin/env python
import numpy as np
from helpers import *
from pcadimreduction import PCA
import parallel as pl
import mprms
import output as oput
from distance import HandleMolCoords as GetCoords
from distance import SetCoords, ScatterCoords
from rdkithelpers import *
#import Coords as crd
import drivers as dr
from glob import glob

RefreshGrid=None
grid=None
fastDecider=True
BinsByVariance=False
nBins=None
DE=True

def Init():
    global nCellDims, grid
    nCellDims = len(nBins)
    if mprms.restart:
        grid=Depickle('grid.P')
    return

#Lay a grid over space in first ndims dimensions, based on boundaries
#of passed data.
#Calls to grid object will return grid assignments
class Grid(object):
    def __init__(self, mols, ndims, nbins):
        coords=self.CoordArray(mols)

        bounds=[ (min(coords[:,i]),max(coords[:,i])) for i in xrange(ndims) ]
        
        self.ndims=ndims
        print "type(nbins):", type(nbins)
        if type(nbins)==int:
            print "nbins==int:", nbins
            print "ndims:", ndims
            print "self.ndims:", self.ndims
            self.nbins=np.array( [nbins]*self.ndims )
        else:
            self.nbins=np.array( nbins )
            assert len(self.nbins)==self.ndims

        self.maxes=np.array( [ max(b) for b in bounds ] )
        self.mins=np.array( [ min(b) for b in bounds ] )
        self.ranges=np.array( [ mx-mn for mx,mn
                                in zip(self.maxes,self.mins) ] )
        self.deltas=np.array( [r/(nbin-1.0)
                               for r, nbin in
                               zip(self.ranges, self.nbins) ] )
        if DE:
            print "self.deltas:", self.deltas
            print "self.ranges:", self.ranges
            print "bounds:", bounds
            print "coords.shape:", coords.shape

        if not hasattr(self,'transform'): self.transform=None


    #Return grid assignments
    def __call__(self,mols):
        coords=self.CoordArray(mols)
        if self.transform is not None:
            pcc = self.transform(coords)[:,:self.ndims]
        else:
            pcc=np.array(coords)

        gridcoords = (pcc[:,:self.ndims]-self.mins)/self.deltas

        if DE:
            print gridcoords
            
        return [ tuple( int(x) for x in gc ) for gc in np.nan_to_num(gridcoords)  ]


    def Draw(self,dimX,dimY, **kwargs):
        from matplotlib.pyplot import hlines,vlines
        minx=self.mins[dimX]
        maxx=self.maxes[dimX]
        miny=self.mins[dimY]
        maxy=self.maxes[dimY]
        deltax=self.deltas[dimX]
        deltay=self.deltas[dimY]

        hlines( [miny+iY*deltay for iY in xrange(self.nbins[dimY]+1)],
                minx, maxx, **kwargs)
        vlines( [minx+iX*deltax for iX in xrange(self.nbins[dimX]+1)],
                miny,maxy,**kwargs)
        

    @staticmethod
    def CoordArray(mols):
        try:
            # I interchanged these two rows
            coords=[SetCoords(m) for m in mols ]
            ndim=len(GetListProp(mols[0], 'coords'))
            #ndim=len(mols[0].GetProp('coords'))
        except AttributeError: coords=mols
        except ValueError:
            # Jos added this. because of:
            # ValueError: GetProp: coords not found.
            coords=mols

        if not issubclass( type(coords), np.ndarray ):
            outarray=np.empty( (len(coords),len(coords[0]) ) )
            outarray[:,:]=coords
        else: outarray=coords
                
        return outarray


#Do a PCA transformation before creating grid
class PCAGrid(Grid):
    def __init__(self, mols, ndims, nbins,
                 pcabasis=None, scaleBins=False):
        coords=self.CoordArray(mols)

        if pcabasis is None:
            pcdecomp = PCA(coords,norm=True)
        else:
            pcdecomp = PCA( self.CoordArray(pcabasis),
                            norm=True )
            
        self.pcdecomp=pcdecomp
        
        pc_coords= pcdecomp.Project(coords)

        self.transform=pcdecomp

        if scaleBins:
            nbinBase=nbins
            evals = pcdecomp.evals
            std_devs = np.sqrt( evals )
            maxdev=std_devs[0]
            nbins=[ int(np.ceil(nbinBase*x/maxdev))
                    for x in std_devs[:ndims] ]
            print 'Grid dimensions: '+'x'.join( map(str,nbins))

        try:
            super(PCAGrid,self).__init__(pc_coords, ndims, nbins[:ndims])
        except TypeError:
            print "nbins:", nbins
            print "ndims:", ndims
            raise
                                     
        
    def Regrid(self, nbins ,scaleBins=False):
        if scaleBins:
            nbinBase=nbins
            evals = self.pcdecomp.evals
            std_devs = np.sqrt( evals )
            maxdev=std_devs[0]
            nbins=[ int(np.ceil(nbinBase*x/maxdev))
                    for x in std_devs[: self.ndims] ]
            print 'Grid dimensions: '+'x'.join( map(str,nbins))

        if type(nbins)==int:
            self.nbins=np.array( [nbins]*self.ndims )
        else:
            self.nbins=np.array( nbins )
            assert len(self.nbins)==self.ndims
            
        self.deltas=np.array( [r/(nbin-1.0)
                               for r,nbin in
                               zip(self.ranges,self.nbins) ] )



#############################################################
#Function to decide which compound stays in;
#Can be replaced
from rdkit.Chem import Descriptors
decider=lambda m:-1.0 * Descriptors.MolWt(m)
lastgrid=None

def GridDiversity(oldmols, newmols,
                  pcabasis=None,molgrid=None):
    global grid, nCellDims, nBins

    if DE:
        print "nCellDims:", nCellDims
        print "nBins:", nBins

    if molgrid is None:
        molgrid = dict()
        mols=oldmols+newmols
    else:
        mols=newmols

    if grid is None:
        print 'Creating new grid'
        grid=PCAGrid( mols,  nCellDims, nBins,
                      pcabasis=pcabasis,
                      scaleBins=BinsByVariance)


    ScatterCoords([m for m in mols if not m.HasProp('gridcoord')])
    ScatterAssign(mols)
    ScatterDecider(mols)

    nNew=0

    oput.StartTimer('GRID PICKS')
    for mol in mols:
        index=mol.GetProp('gridcoord')
        if molgrid.has_key( index ):
            if molgrid[index].GetDoubleProp('decider') < mol.GetDoubleProp('decider'):
                molgrid[index]=mol
                nNew+=1
        else:
            molgrid[index]=mol
            nNew+=1

    print len(molgrid),'/',np.product(grid.nbins),'occupied cells (',\
          nNew,'new)'
    oput.EndTimer('GRID PICKS')
            
    return len(molgrid) , molgrid.values(),molgrid

    

#######################################################################
#Grid diversity that delays filtering until AFTER looking for diversity
#def GridDiversity_JITFilter(oldmols,newmols,Filter,Geom,molgrid=None):
def GridDiversity_JITFilter(oldmols, newmols, Filter=True, Geom=False, molgrid=None):

    if DE:
        print "oldmols:", oldmols

    if not (Filter or Geom):
        return GridDiversity( oldmols, newmols )

    if Filter: newmols=dr.DriveFilters( newmols, Filter, False)
    ScatterCoords([m for m in newmols if not m.HasProp('gridcoord')])
    ScatterAssign(newmols)
    ScatterDecider(newmols)

    #Get old assignments if not passed
    oput.StartTimer('GRID PICKS')
    if not molgrid:
        molgrid={}
        for mol in oldmols:
            index=mol.GetProp('gridcoord')
            if molgrid.has_key( index ):
                if molgrid[index].GetDoubleProp('decider') < mol.GetDoubleProp('decider'):
                    molgrid[index]=mol
            else:
                molgrid[index]=mol

    nOld=len(molgrid)

    #Screen new molecules for novelty
    toscreen=[]
    for mol in newmols:
        index=mol.GetProp('gridcoord')
        try:
            if (not molgrid.has_key(index)) or \
                molgrid[index].GetDoubleProp('decider')<mol.GetDoubleProp('decider'):
                toscreen.append(mol)
        except ValueError:
            print "No Decider keyword for:", Chem.MolToSmiles(molgrid[index])
            #print molgrid[index].GetDoubleProp('decider')
            pass

    oput.EndTimer('GRID PICKS')

    print 'Novel mutants:',len(toscreen),'/',len(newmols)

    #Run filters on novel molecules only
    if Geom:
        goodmols=dr.DriveFilters(toscreen,Filter,Geom)
        print 'Molecules passing filters:',len(goodmols)
    else:
        goodmols=toscreen

    #Check for some problems
    badgoodmols=[m for m in goodmols if not m.HasProp('gridcoord')]
    if len(badgoodmols)>0:
        DumpLowMemMols(badgoodmols,'failchange.pjar.gz',True)
        print 'Molecule was changed after filtering ...'
        ScatterCoords(badgoodmols)
        ScatterAssign(badgoodmols)
        ScatterDecider(badgoodmols)


    #Assign filtered novel molecules
    oput.StartTimer('GRID PICKS')
    nReplace=0
    for mol in goodmols:
        index=mol.GetProp('gridcoord')
        if molgrid.has_key( index ):
            if molgrid[index].GetDoubleProp('decider') < mol.GetDoubleProp('decider'):
                molgrid[index]=mol
                nReplace+=1
        else:
            molgrid[index]=mol
    oput.EndTimer('GRID PICKS')

    print len(molgrid),'/',np.product(grid.nbins),'occupied cells (',\
          len(molgrid)-nOld,'new, '+str(nReplace)+' replaced.)'
        
    return len(molgrid),molgrid.values(),molgrid
    



ChunkSize=200000
def ScatterAssign(mols):
    toassign=[m for m in mols if not m.HasProp('gridcoord')]
    if len(toassign)==0: return
    
    oput.StartTimer('GRID ASSIGN')
    if len(toassign)<ChunkSize or not pl.mpi:
        assignments=grid(toassign)
        for assign,mol in zip(assignments,toassign):
            print "type(assign):", type(assign)
            print "assign:", assign
            print [ type(item) for item in assign ]
            newassign=tuple(map(float,assign))
            SetListProp(mol, 'gridcoord', newassign)
            

    else:
        print 'Scattering grid assignments ...'
        scatter=[ (grid,
                   [m.GetProp('coords') for m in toassign[i:i+ChunkSize] ] )
                  for i in xrange(0,len(toassign),ChunkSize ) ]
        if pl.verbose:
            print len(scatter)
        
        
        pl.MyTask.SetFunction(MPIGridAssign)
        assigns=pl.MyTask.RunMPI(scatter)
        for mol, gridcoord in zip(toassign,(a for x in assigns for a in x)):
            mol.SetProp( 'gridcoord' , gridcoord )
        
    
    oput.EndTimer('GRID ASSIGN')
    

#@pl.MPIScatter
def MPIGridAssign(args):
    grid,coords=args

    assigns=grid(coords)

    return assigns


moltemp=Chem.Mol()
def ScatterDecider(mols):
    toCompute=[m for m in mols if not m.HasProp('decider')]
    if len(toCompute)==0: return
    else:
        pass
        #print "len(toCompute):", len(toCompute)
        #print 'C1=CC(=O)C(=O)C=C1CC=O' in [ oe.OECreateCanSmiString(m) for m in mols ]

    oput.StartTimer('BIAS FUNCTION')
    if not pl.mpi or fastDecider:
        for m in toCompute:
            print m.SetDoubleProp('decider',decider(m))
    else:
        print 'Scattering bias function ...'
        sendmols = [pl.SendMol(m) for m in toCompute]
        pl.MyTask.SetFunction(MPIDecider)
        vals=pl.MyTask.RunMPI( sendmols )
        for v,m in zip(vals,toCompute):
            m.SetDoubleProp('decider' , v )
    oput.EndTimer('BIAS FUNCTION')


#@pl.MPIScatter
def MPIDecider(m):
    return decider( pl.RecvMol(m, oemol=True) )
        
        
WatchFolder=None
WatchPrefix='it'
NoReadFiles=set()
mymaxit=-1
def ReadWatchFolder(mylib,wrotepool):
    global mymaxit

    oput.StartTimer('READ')
    if pl.mpi:
        answer=ScatterReadWatchFolder(mylib,wrotepool)
        oput.EndTimer('READ')
        return answer

    mysmi=set( m.GetProp('isosmi') for m in mylib )
    
    print "Reading additional molecules: ",
    for file in glob( WatchFolder+'*.oeb.gz'):

        if basename(file)[:len(WatchPrefix)]==WatchPrefix:
            mymaxit = extractnum(basename(file))

        if basename(file,False) in NoReadFiles: continue

        nNew=0
        for nmol,newmol in enumerate(GetLowMemMols(file)):
            smi=newmol.GetProp('isosmi')
            if not smi in mysmi:
                nNew+=1
                mysmi.add( smi )
                wrotepool.add( smi )
                StripData( newmol )
                mylib.append( newmol )
                
        NoReadFiles.add( basename(file,False) )
        print basename(file,True)+' ('+str(nNew)+'/'+str(nmol+1)+')',

    print 'done.'
    oput.EndTimer('READ')
    return mymaxit


#MPI version of above
def ScatterReadWatchFolder(mylib,wrotepool):
    global mymaxit
    import random

    mysmi=set( m.GetProp('isosmi') for m in mylib )
    print 'Reading additional molecules'

    #Get new files to read in
    filesToRead=list( set( glob(WatchFolder+'*.oeb.gz') )
                      - NoReadFiles )
    NoReadFiles.update( filesToRead )

    #Make sure one node isn't stuck with all the big files
    random.shuffle(filesToRead)

                      
    if len(filesToRead)==0: return mymaxit

    pl.MyTask.SetFunction(MPIReadFiles)
    readMols=pl.MyTask.RunMPI( filesToRead )

    newmols=[]
    for fname,mlist in zip(filesToRead,readMols):

        if basename(fname)[:len(WatchPrefix)]==WatchPrefix:
            mymaxit = extractnum(basename(fname))
            
        nNew=0
        for nmol,newmol in enumerate(mlist):
            smi=newmol.GetProp('isosmi')
            if not smi in mysmi:
                nNew+=1
                mysmi.add( smi )
                wrotepool.add( smi )
                StripData( newmol )
                mylib.append( newmol )
        print basename(fname,True)+' ('+str(nNew)+'/'+str(nmol+1)+')',
                       
    print 'done.'
    return mymaxit
        

    


#@pl.MPIScatter
def MPIReadFiles(fname):
    result=[]
    for m in GetLowMemMolIter(fname):
        StripData(m)
        result.append(m)
    print 'File: '+fname,len(result)
    return result
    

def WriteWatchFolder(mylib,wrotepool,itnum):
    oput.StartTimer('WRITE')
    filename=WatchPrefix+str(itnum) + '.oeb.gz'
    writemols=list( LMMScreener(mylib,wrotepool) )
    DumpMols( writemols , filename)
    for mol in writemols: StripData(mol)
    os.system('mv '+filename+' '+WatchFolder+filename)
    NoReadFiles.add(WatchFolder+filename)
    oput.EndTimer('WRITE')



#For really low memory bookkeeping, this is all we need
keepers=set(('filtered','failedfilter','hasstructure',
             'isosmi','gridcoord','decider','Objective',
             'selected'))
def StripData(mol):
    delDatas=set(mol.GetProp().keys())-keepers
    for dat in delDatas: mol.DeleteData(dat)
        
            
    
def NewListsFromAssignments(mygrid,coords,**kwargs):
    if kwargs.has_key('decider'):
        decider=kwargs['decider']
    else:
        decider=None

    outnames=kwargs.keys()
    inlists=[kwargs[key] for key in outnames]

    outnames.append('coords')
    inlists.append(coords)

    outgrids={ name:{} for name in outnames }

    assigns=mygrid(coords)

    outnames.insert(0,'gridcoords')
    inlists.insert(0,assigns)

#    for vals in zip( *inlists ):
#        if

    
