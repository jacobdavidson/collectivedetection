import numpy as np
import scipy.spatial

def fishmodel(pos,le,re,tail):
    # pos = head, le=left eye, re = right eye, tail = tail
    
    points = np.array([pos.T,re.T,tail.T,le.T,pos.T])
      
    # make a 4 sided figure
    numparticles=len(pos)
    fishsegments=np.zeros([numparticles,4,2,2])

    for i in range(4):
        fishsegments[:,i,0]=points[i].T
        fishsegments[:,i,1]=points[np.mod(i+1,4)].T
    return fishsegments  ## segments 0 and 1 are right eye, 2 and 3 are left eye

def makesegmentsinglestep(pos,head,fishlength):
    numparticles=len(pos)
    fishsegments=np.zeros([numparticles,2,2])
    for fish in range(numparticles):
        fishsegments[fish]=[[pos[fish,0]+0*fishlength/2*np.cos(head[fish]),pos[fish,1]+0*fishlength/2*np.sin(head[fish])],
                    [pos[fish,0]-fishlength*np.cos(head[fish]),pos[fish,1]-fishlength*np.sin(head[fish])]]
    return fishsegments

def getpointpositions(xy,hxhy,dist):
    x, y = xy
    hx, hy = hxhy
    front = [x+hx*dist, y+hy*dist]
    back = [x-hx*dist, y-hy*dist]
    left = [x-hy*dist, y+hx*dist] # not sure which is left and which is right
    right = [x+hy*dist, y-hx*dist]
    return np.array([front,back,left,right])

def getpointpositions_many(xy,hxhy,dist,numpoints=32):
    x, y = xy
    hx, hy = hxhy
    dirangle = np.arctan2(hy,hx)
    angles = np.linspace(0,2*np.pi,numpoints+1)[0:-1] + dirangle
    positions = np.array([x+dist*np.cos(angles),y+dist*np.sin(angles)]).T
    return positions

def fishlines_to_segments(xp,yp):
    segments=np.zeros([len(xp)-1,2,2])
    for k in range(len(xp)-1):
        segments[k] = [[xp[k],yp[k]],[xp[k+1],yp[k+1]]]
    return segments

# get segments                  
def getbonesegsinglestep(pos,midc,skip=1):
    segstep=[]
    numfish, numbones, _ = midc.shape 
    toselect = np.arange(0,numbones,skip)
    if ~(toselect[-1]==numbones-1):  # make sure the last point in the list
        toselect = np.append(toselect,numbones-1)
    for i in range(numfish):
        xb = pos[i,0]+midc[i,toselect,0]
        yb = pos[i,1]+midc[i,toselect,1]
        bseg=fishlines_to_segments(xb,yb)       
        segstep.append(bseg)
    return np.array(segstep).reshape(-1,2,2)

# raycasting-related functions
def fixanglerange(angles): # Puts all angles into a range of [-Pi,Pi] 
    return np.arctan2(np.sin(angles),np.cos(angles))


def getoneray_intersecttest_nonopt(viewpoint,segments,rayangles):  # this is a slower way of doing the same thing
    numsegments=len(segments)

    intersect = np.zeros(len(rayangles))
    for seg in segments:
        if np.sum(intersect)<len(intersect):  # skip if all already are found to intersect
            diffs=np.array([seg[0]-viewpoint,seg[1]-viewpoint])
            thetas=list(map(lambda x: np.arctan2(x[1],x[0]),diffs))        
            thetadiff=fixanglerange(thetas[1]-thetas[0])
            # always use a positive angle difference:  reverse points if angle is negative
            if thetadiff<0:
                #distances=np.flip(distances,0)
                thetas=np.flip(thetas,0)
                thetadiff=-thetadiff
            # thetadiff will now always be nonnegative.  make sure its nonzero, to prevent errors, in case its zero
            if thetadiff>0:
                for rnum in range(len(rayangles)):
                    rayangle=rayangles[rnum]
                    thetastart = thetas[0].copy()
                    # correct for the wrap-around case
                    if thetas[0]-rayangle>np.pi: 
                        thetastart=thetastart-2*np.pi
                    if thetas[0]-rayangle<-np.pi: 
                        thetastart=thetastart+2*np.pi

                    # now, check if its within this range
                    if (rayangle>=thetastart) & (rayangle<=thetastart+thetadiff):
                        intersect[rnum] = 1
    return np.logical_not(intersect.astype(bool))

def getoneray_intersecttest(viewpoint,segments,rayangles):
    numsegments=len(segments)

    intersect = np.zeros(len(rayangles))
    
    alldiffs = segments-viewpoint
    allthetas = np.arctan2(alldiffs[:,:,1],alldiffs[:,:,0])
            
    for thetas in allthetas:
        thetadiff=fixanglerange(thetas[1]-thetas[0])
        # always use a positive angle difference:  reverse points if angle is negative
        if thetadiff<0:
            #distances=np.flip(distances,0)
            thetas=np.flip(thetas,0)
            thetadiff=-thetadiff
        # thetadiff will now always be nonnegative.  make sure its nonzero, to prevent errors, in case its zero
        if thetadiff>0:
            for rnum in range(len(rayangles)):
                if intersect[rnum]==0:
                    rayangle=rayangles[rnum]
                    thetastart = thetas[0].copy()
                    # correct for the wrap-around case
                    if thetas[0]-rayangle>np.pi: 
                        thetastart=thetastart-2*np.pi
                    if thetas[0]-rayangle<-np.pi: 
                        thetastart=thetastart+2*np.pi

                    # now, check if its within this range
                    if (rayangle>=thetastart) & (rayangle<=thetastart+thetadiff):
                        intersect[rnum] = 1
    return np.logical_not(intersect.astype(bool))


def getoneray_intersecttest_probabilistic(viewpoint,segments,rayangles,randvalues,blockingprobability=1):
    # randvalues should have the shape [len(segments),len(rayangles)].  But pre-generate for left and right eye consistency
    numsegments=len(segments)

    intersect = np.zeros(len(rayangles))
    
    alldiffs = segments-viewpoint
    allthetas = np.arctan2(alldiffs[:,:,1],alldiffs[:,:,0])
            
    for t in range(len(allthetas)):
        thetas = allthetas[t]
        thetadiff=fixanglerange(thetas[1]-thetas[0])
        # always use a positive angle difference:  reverse points if angle is negative
        if thetadiff<0:
            #distances=np.flip(distances,0)
            thetas=np.flip(thetas,0)
            thetadiff=-thetadiff
        # thetadiff will now always be nonnegative.  make sure its nonzero, to prevent errors, in case its zero
        if thetadiff>0:
            for rnum in range(len(rayangles)):
                if intersect[rnum]==0:
                    rayangle=rayangles[rnum]
                    thetastart = thetas[0].copy()
                    # correct for the wrap-around case
                    if thetas[0]-rayangle>np.pi: 
                        thetastart=thetastart-2*np.pi
                    if thetas[0]-rayangle<-np.pi: 
                        thetastart=thetastart+2*np.pi

                    # now, check if its within this range
                    if (rayangle>=thetastart) & (rayangle<=thetastart+thetadiff):
                        if randvalues[t,rnum]<blockingprobability:
                            intersect[rnum] = 1
    return np.logical_not(intersect.astype(bool))

def getvoronai(positions):
    vnet=np.zeros((len(positions),len(positions)))
    tri = scipy.spatial.Delaunay(positions)
    indices,indptr=tri.vertex_neighbor_vertices                  
    for i in range(len(positions)):
        connections=np.sort(indptr[indices[i]:indices[i+1]])  # see here :  https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.spatial.Delaunay.vertex_neighbor_vertices.html#scipy.spatial.Delaunay.vertex_neighbor_vertices
        vnet[i,connections]=1
    return vnet