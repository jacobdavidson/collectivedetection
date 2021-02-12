import numpy as np

# raycasting functions
def vectorlength(x):
    return np.sqrt(x.dot(x))
def fixanglerange(angles): # Puts all angles into a range of [-Pi/2,Pi/2] 
    return np.arctan2(np.sin(angles),np.cos(angles))


def getrayrangeindices(thetas,numrays):
    dray = 2*np.pi/numrays
    return np.array([np.mod(np.ceil((thetas[0]+np.pi)/dray),numrays),
                      np.floor((thetas[1]+np.pi)/dray)]
                     ,dtype=np.int)

def checkrange(rayrange,th):
    newrange=rayrange
    if (rayrange[0]<-np.pi) | (rayrange[1]>np.pi) :
        newrange=np.sort([fixanglerange(r) for r in rayrange])
        return np.logical_not((th>=newrange[0]) & (th<=newrange[1]))
    else:
        return (th>=newrange[0]) & (th<=newrange[1])

def getrays(viewpoint,segments,numrays,rayrange=[-np.pi,np.pi],ignore=[]):
    numsegments=len(segments);
    rayvalues = 2*np.pi*np.arange(0,numrays)/numrays - np.pi
    dray = rayvalues[1] - rayvalues[0]
    raysingle = -1*np.ones((numrays, 2))  # for each ray, stores [lnum of intersect, distance]
    
    for lnum in range(numsegments):
        if not (lnum in ignore):
            diffs=np.array([segments[lnum,0]-viewpoint,segments[lnum,1]-viewpoint])
            # don't actually need this!
            #distances=list(map(vectorlength,diffs))
            thetas=list(map(lambda x: np.arctan2(x[1],x[0]),diffs))
            thetadiff=fixanglerange(thetas[1]-thetas[0])
            # always use a positive angle difference:  reverse points if angle is negative
            if thetadiff<0:
                #distances=np.flip(distances,0)
                thetas=np.flip(thetas,0)
                thetadiff=-thetadiff
            # thetadiff will now always be nonnegative.  make sure its nonzero, to prevent errors
            if thetadiff>0:
                rangeindices=getrayrangeindices(thetas,numrays)
                # numindices is the number of indices that it overlaps with - only check these, not the others
                numindices=np.round(fixanglerange((rayvalues[rangeindices[1]]-rayvalues[rangeindices[0]]))/dray)+1
                numindices=numindices.astype(int)
                
                ### JD, on Oct 29, 2018:  I am not sure what this does here, with 'thetastart'.  The code doesn't seem to use this
                thetastart=thetas[0]
                # correct for the wrap-around case
                if thetas[0]-rayvalues[rangeindices[0]]>np.pi: 
                    thetastart=thetastart-2*np.pi
                if thetas[0]-rayvalues[rangeindices[0]]<-np.pi: 
                    thetastart=thetastart+2*np.pi
                # another wrap-around correction
                if numindices>numrays/2:
                    numindices=0

                # the distance is calculated by solving for the intersection point.  that is how get this formula
                # these give the same thing
                rdistcommon=np.abs((-diffs[0,1] + diffs[1,1])*viewpoint[0] + (diffs[0,0]-diffs[1,0])*viewpoint[1] +  segments[lnum,0,1]*segments[lnum,1,0]-segments[lnum,0,0]*segments[lnum,1,1] )
                # do the simplest way:  check all entries (could be faster to check first for conflicts, and if not, replace all, but this is simpler)
                for rnum in range(numindices):
                    rindex=rangeindices[0]+np.sign(thetadiff)*(rnum)
                    rindex=np.mod(rindex,numrays).astype(int)
                    if checkrange(rayrange,rayvalues[rindex]):
                        #calculate distance to intersection  
                        # these give the same thing, use either one
                        rdist=rdistcommon/np.abs((diffs[0,1]-diffs[1,1])*np.cos(rayvalues[rindex]) - (diffs[0,0] - diffs[1,0])*np.sin(rayvalues[rindex]))          
                        #rdist=rdistcommon/np.abs(np.cos(rayvalues[rindex])*(diffs[0,1]-diffs[1,1] + (-diffs[0,0]+diffs[1,0])*np.tan(rayvalues[rindex]) ))
                        if raysingle[rindex,0]==-1:  # then fill in, there is nothing there yet
                            raysingle[rindex,0]=lnum
                            raysingle[rindex,1]=rdist
                        elif rdist<raysingle[rindex,1]:  # else if there is an entry already there, replace it if the new ray is closer
                            raysingle[rindex,0]=lnum
                            raysingle[rindex,1]=rdist   
    return raysingle

