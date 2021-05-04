#!/usr/bin/gawk -f
#
# units:  cm, radian
#         cm, degree (on command line)
#
# The coordinate system: z-beam direction, x-up, y-right (looking in beam direction).
# The original crystal location: centered at origin, #
#                                smallest face in z-direction
#                                largest face in x-direction
# coordinates are named as follows
# name |  x   |  y   |  z   |  comment
# +++: |+sx1/2|+sy1/2| sz/2 | upper-right corner of downstream face
# ++-: |+sx2/2|+sy2/2|-sz/2 | upper-right corner of upstream face
# +0-: |+sx2/2|  0   |-sz/2 | middle of right upstream edge
# c    |                      center
# xd                          initial (or intrinsic) x-direction
# yd                          ...

#
# 2008-05-16  HS  renamed to 'det_place.gawk'  (detector placement) 
# 2008-05-20  HS  major rewrite, incluing trapezoidal detectors
# 2008-09-15  HS  added option to place detectors in groups of 2,3,4... with sides in theta-direction joined parallel touching
# 2008-10-10  HS  added option to specify wall thickness (-w)

function setupTypeSpec(t) {
# detector specs
# type 1: new DALI
    t = 1;
    typeSpec[t,"sx1"] = 4.5;   # size in x, y, z of a trapezoid
    typeSpec[t,"sx2"] = 4.5;   # sx1 is the size at z=+sz/2
    typeSpec[t,"sy1"] = 8.3;   # sx2 is the size at z=-sz/2
    typeSpec[t,"sy2"] = 8.3;   # same for y
    typeSpec[t,"sz"] = 18;
# type 2: new DALI
    t = 2;
    typeSpec[t,"sx1"] = 5.01;   # size in x, y, z
    typeSpec[t,"sx2"] = 5.01;   
    typeSpec[t,"sy1"] = 8.52;
    typeSpec[t,"sy2"] = 8.52;
    typeSpec[t,"sz"] = 17.64;
# type 3: old DALI
    t = 3;
    typeSpec[t,"sx1"] = 6.6;   # size in x, y, z
    typeSpec[t,"sx2"] = 6.6;   
    typeSpec[t,"sy1"] = 6.6;
    typeSpec[t,"sy2"] = 6.6;
    typeSpec[t,"sz"] = 12.6;
# type 4: possible LaBr or LaCl size
    t = 4;
    typeSpec[t,"sx1"] = 1.5;   # size in x, y, z  (suggested LaBr size)
    typeSpec[t,"sx2"] = 1.5;   
    typeSpec[t,"sy1"] = 4;
    typeSpec[t,"sy2"] = 4;
    typeSpec[t,"sz"] = 6;
# type 5: possible LaBr or LaCl size
    t = 5;
    typeSpec[t,"sx1"] = 1.5;   # size in x, y, z  (somewhat longer than type 4)
    typeSpec[t,"sx2"] = 1.5;   
    typeSpec[t,"sy1"] = 4;
    typeSpec[t,"sy2"] = 4;
    typeSpec[t,"sz"] = 8;
# type 14: possible LaBr or LaCl size
    t = 14;
    typeSpec[t,"sx1"] = 2;   # size in x, y, z  (suggested LaBr size)
    typeSpec[t,"sx2"] = 1.5;   
    typeSpec[t,"sy1"] = 5;
    typeSpec[t,"sy2"] = 4;
    typeSpec[t,"sz"] = 6;
# type 15: possible LaBr or LaCl size
    t = 15;
    typeSpec[t,"sx1"] = 2;   # size in x, y, z  (somewhat longer than type 4)
    typeSpec[t,"sx2"] = 1.5;   
    typeSpec[t,"sy1"] = 5;
    typeSpec[t,"sy2"] = 4;
    typeSpec[t,"sz"] = 8;
# type 24: possible LaBr or LaCl size
    t = 24;
    typeSpec[t,"sx1"] = 1.5;   # size in x, y, z  (suggested LaBr size)
    typeSpec[t,"sx2"] = 1.5;   
    typeSpec[t,"sy1"] = 5;
    typeSpec[t,"sy2"] = 4;
    typeSpec[t,"sz"] = 6;
# type 25: possible LaBr or LaCl size
    t = 25;
    typeSpec[t,"sx1"] = 1.5;   # size in x, y, z  (somewhat longer than type 4)
    typeSpec[t,"sx2"] = 1.5;   
    typeSpec[t,"sy1"] = 5;
    typeSpec[t,"sy2"] = 4;
    typeSpec[t,"sz"] = 8;
# type 6: possible CsI size
    t = 6;
    typeSpec[t,"sx1"] = 2.5;   # size in x, y, z
    typeSpec[t,"sx2"] = 2.5;   
    typeSpec[t,"sy1"] = 5;
    typeSpec[t,"sy2"] = 5;
    typeSpec[t,"sz"] = 16;
# type 501: trapezoidal for testing
    t = 501;
    typeSpec[t,"sx1"] = 4;   # size in x, y, z
    typeSpec[t,"sx2"] = 2;   
    typeSpec[t,"sy1"] = 8;
    typeSpec[t,"sy2"] = 6;
    typeSpec[t,"sz"] = 20;

# type 555: rather small size for testing
    t = 555;
    typeSpec[t,"sx1"] = 1.0;   # size in x, y, z
    typeSpec[t,"sx2"] = 1.0;   
    typeSpec[t,"sy1"] = 1.0;
    typeSpec[t,"sy2"] = 1.0;
    typeSpec[t,"sz"] = .1;
}

#
# tan(x)
function tan(x) {
    return cos(x) != 0 ? sin(x)/cos(x) : sin(x)/0;
}

#
# reset everything for detector i
#
function resetDet(i) {
        det[i,"r"]  = 0;     # ring ID
        det[i,"d"]  = 0;     # distance
        det[i,"th"] = 0;     # theta
        det[i,"ph"] = 0;     # phi
        det[i,"u"]  = 0;     # used (flag if used or not)
        resetPos(i);
}

#
# reset crystal position and orientation
#
function resetPos(i, \
                  c, w, x, xx, tok, f, sx, sy, sz, ii) { 
    for (w in direct) {
        for (c in coord) {
            det[i,c,w]= 0;
        }
    }
    # initialize direction
    det[i,"xd","x"] = 1;   # intrinsic x-direction in the laboratory (i.e. before placing the detectors)
    det[i,"yd","y"] = 1;   # intrinsic y-direction in the laboratory
    det[i,"zd","z"] = 1;   # intrinsic z-direction in the laboratory
# properly populate all +++ coordinates
    for (c in coord) {
        if (c !~ /[+0-][+0-][+0-]/) continue;
        for (ii=1; ii<4; ii++) {
            tok = substr(c,ii,1);
            f[ii] = tok 1;
            if (tok == "0") f[ii]=0;
            f[ii] *= 1.0
        }
        sz = det[i,"sz"]/2;
        if ( 0 ) {
        } else if (f[3] > 0 ) {
            sx = det[i,"sx1"]/2;
            sy = det[i,"sy1"]/2;
        } else if (f[3] < 0 ) {
            sx = det[i,"sx2"]/2;
            sy = det[i,"sy2"]/2;
        } else if (f[3] == 0 ) {
            sx = det[i,"sx2"]/2 + det[i,"sx1"]/2;
            sy = det[i,"sy2"]/2 + det[i,"sy1"]/2;
            sx *= .5;
            sy *= .5;
        } else {
            print "unknown value of f[3]";
            exit 1;
        }
        det[i,c,"x"] = sx*f[1];
        det[i,c,"y"] = sy*f[2];
        det[i,c,"z"] = sz*f[3];
    }
}

#
#
#
function showDet (i, \
                  x, xx) {
    for (x in det) {
        split(x, xx, "\034");
        if (xx[1] == i) {
            printf "%4s %6s %2s %s\n",xx[1], xx[2], xx[3], det[x]
        }
    }
}

# "arguments" on second line are local variables
function thetaPhiCover (i, ret, \
                        c, t, p, th, tl, tc, ph, pl, pc, xd, yd, zd, firsttime) {
# In this routine the lowest and highest theta and phi values are determined
    firsttime = 1;
    for (c in coord) {
        if (c !~ /[+0-][+0-][+0-]/) continue;
        xd = det[i, c, "x"];
        yd = det[i, c, "y"];
        zd = det[i, c, "z"];
#        print "xd yd zd",xd,yd,zd
        t = atan2(sqrt(xd^2 + yd^2),zd);
        p = atan2(yd,xd);
#        print " i, t, p  " c, t, p
        # save theta and phi of each point
        det[i, c, "theta"] = t;
        det[i, c, "phi"]   = p;
        if (firsttime) {
            th = t;
            tl = t;
            ph = p;
            pl = p;
            firsttime=0;
            continue;
        }
        if ( t > th ) th = t;
        if ( t < tl ) tl = t;
        if ( p > ph ) ph = p;
        if ( p < pl ) pl = p;
    }
    tc = th - tl;  # theta cover
    pc = ph - pl;  # phi cover
#    print "th tl tc ", th/DEG, tl/DEG, tc/DEG
#    print "ph pl pc ", ph/DEG, pl/DEG, pc/DEG

# shape in xz plane (middle of detector)
    sxz = "";
    xd = det[i, "+0+", "x"];
    zd = det[i, "+0+", "z"];
    sxz = sxz zd " " xd "\n";
    xd = det[i, "+0-", "x"];
    zd = det[i, "+0-", "z"];
    sxz = sxz zd " " xd "\n";
    xd = det[i, "-0-", "x"];
    zd = det[i, "-0-", "z"];
    sxz = sxz zd " " xd "\n";
    xd = det[i, "-0+", "x"];
    zd = det[i, "-0+", "z"];
    sxz = sxz zd " " xd "\n";
    xd = det[i, "+0+", "x"];
    zd = det[i, "+0+", "z"];
    sxz = sxz zd " " xd "\n";

# shape in xy plane (looking into the beam x = y-axis; -y = x-axis)
    sxy = "";
# back face
    xd = det[i, "+++","x"];
    yd = det[i, "+++","y"];
    sxy = sxy yd " " xd "\n";
    xd = det[i, "+-+","x"];
    yd = det[i, "+-+","y"];
    sxy = sxy yd " " xd "\n";
    xd = det[i, "--+","x"];
    yd = det[i, "--+","y"];
    sxy = sxy yd " " xd "\n";
    xd = det[i, "-++","x"];
    yd = det[i, "-++","y"];
    sxy = sxy yd " " xd "\n";
    xd = det[i, "+++","x"];
    yd = det[i, "+++","y"];
    sxy = sxy yd " " xd "\n";

# front face
    sxy = sxy "\n";
    xd = det[i,"++-","x"]
    yd = det[i,"++-","y"]
    sxy = sxy yd " " xd "\n";
    xd = det[i,"+--","x"]
    yd = det[i,"+--","y"]
    sxy = sxy yd " " xd "\n";
    xd = det[i,"---","x"]
    yd = det[i,"---","y"]
    sxy = sxy yd " " xd "\n";
    xd = det[i,"-+-","x"]
    yd = det[i,"-+-","y"]
    sxy = sxy yd " " xd "\n";
    xd = det[i,"++-","x"]
    yd = det[i,"++-","y"]
    sxy = sxy yd " " xd "\n";


# lower face
    sxy = sxy "\n";
    xd = det[i,"-++","x"];
    yd = det[i,"-++","y"];
    sxy = sxy yd " " xd "\n";
    xd = det[i,"-+-","x"];
    yd = det[i,"-+-","y"];
    sxy = sxy yd " " xd "\n";
    xd = det[i,"---","x"];
    yd = det[i,"---","y"];
    sxy = sxy yd " " xd "\n";
    xd = det[i,"--+","x"];
    yd = det[i,"--+","y"];
    sxy = sxy yd " " xd "\n";
    xd = det[i,"-++","x"];
    yd = det[i,"-++","y"];
    sxy = sxy yd " " xd "\n";

# upper face
    sxy = sxy "\n";
    xd = det[i,"+++","x"];
    yd = det[i,"+++","y"];
    sxy = sxy yd " " xd "\n";
    xd = det[i,"++-","x"];
    yd = det[i,"++-","y"];
    sxy = sxy yd " " xd "\n";
    xd = det[i,"+--","x"];
    yd = det[i,"+--","y"];
    sxy = sxy yd " " xd "\n";
    xd = det[i,"+-+","x"];
    yd = det[i,"+-+","y"];
    sxy = sxy yd " " xd "\n";
    xd = det[i,"+++","x"];
    yd = det[i,"+++","y"];
    sxy = sxy yd " " xd "\n";

#  print sxy

    det[i,"tp-tl"] = tl;  # theta low
    det[i,"tp-th"] = th;  # theta high
    det[i,"tp-tc"] = tc;  # theta cover
    det[i,"tp-pl"] = pl;  # phi low
    det[i,"tp-ph"] = ph;  # phi high
    det[i,"tp-pc"] = pc;  # phi cover
    det[i,"tp-shape xz"] = sxz;
    det[i,"tp-shape xy"] = sxy;

# return values depending on second argument
    if ( ret == "" ) print "cover", tc/DEG, pc/DEG;
    if ( ret == "" ) return 0;
#    if ( ret == "tl" ) return tl;  # theta low
#    if ( ret == "th" ) return th;  # theta high
#    if ( ret == "tc" ) return tc;  # theta center
#    if ( ret == "pl" ) return pl;  # phi low
#    if ( ret == "ph" ) return ph;  # phi high
#    if ( ret == "pc" ) return pc;  # phi center
#    if ( ret == "shape xz" ) return sxz;
#    if ( ret == "shape xy" ) return sxy;
    return 0;
}

#
# initialize detector i to be type t
#
function initDet(i, t) {
    if ( typeSpec[t,"sx1"] == "" || typeSpec[t,"sy1"] == "" || typeSpec[t,"sz"] == "" || \
         typeSpec[t,"sx2"] == "" || typeSpec[t,"sy2"] == "" ){
        print "type " t " not found";
        exit 1;
    }
    det[i,"t"]   = t;                   # type
    det[i,"sx1"] = typeSpec[t,"sx1"];   # size 
    det[i,"sx2"] = typeSpec[t,"sx2"];  
    det[i,"sy1"] = typeSpec[t,"sy1"];
    det[i,"sy2"] = typeSpec[t,"sy2"];
    det[i,"sz"]  = typeSpec[t,"sz"];
}

#
#
# local variables on second line
# rotation around d-axis of detector number n
function rotDet (i, d, v, \
                 c, s, cc, raxis)  {   
    c = cos(v);
    s = sin(v);
    # now do the rotation
    for (cc in coord) {
        for (raxis=1; raxis<4; raxis++) {
            if ( 0 ) {
            } else if (raxis == 1) {                 
                d1="z"; d2="x"; d3="y";   # rot around z-axis
            } else if (raxis == 2) {                 
                d1="x"; d2="y"; d3="z";   # rot around x-axis
            } else if (raxis == 3) {                 
                d1="y"; d2="z"; d3="x";   # rot around y-axis
            } else {
                print "rotDet: unknown raxis"
            }
            if (d != d1)     continue;
            ox = det[i, cc, d2];
            oy = det[i, cc, d3];
            det[i, cc, d2] = ox * c - oy * s;
            det[i, cc, d3] = ox * s + oy * c;
        }
    }
}

function transDet (i, d, v, \
                   c) {   # translation of detector n in direction d by amount v
    # this affects all coordinates, except the direction giving ones
    for (c in coord) {
        if (c ~ /d$/) continue  # if the second index ends with 'd' then it's a direction
        det[i, c, d] += v;
    }
}

function placeDet (i, d, th, ph, x0) {  # place detector i at distance d and angles theta and phi
    transDet(i,"z",d);
    transDet(i,"x",x0);
    rotDet(i,"y",th);
    rotDet(i,"z",ph);
}

#
#
# setup of detector array det[]
function detSetup(i, t, n) {
    for (i = 1; i <= numType; i++) {
        t = type[i,"t"];
        n = type[i,"n"];
        for (j = 1; j <= n; j++) {  # number of detectors of type t
            numDet++;
            initDet(numDet,t);
        }
    }
}

function cvsOutPut() {
    # output of type, positions (x,y,z) and angles (theta, phi) as CSV
    posOut = "t.dat"  # output file for positions
    for (i=1; i<=numDet; i++) {
        if ( ! det[i,"u"] ) continue;
        if (i!=1) printf "," >posOut;
        printf "%d",det[i,"t"] >posOut;
    }
    printf "\n" >posOut;
    
    for (i=1; i<=numDet; i++) {
        if ( ! det[i,"u"] ) continue;
        if (i!=1) printf "," >posOut;
        printf "%.3f",det[i,"x"] >posOut;
    }
    printf "\n" >posOut;
    
    for (i=1; i<=numDet; i++) {
        if ( ! det[i,"u"] ) continue;
        if (i!=1) printf "," >posOut;
        printf "%.3f",det[i,"y"] >posOut;
    }
    printf "\n" >posOut;
    
    for (i=1; i<=numDet; i++) {
        if ( ! det[i,"u"] ) continue;
        if (i!=1) printf "," >posOut;
        printf "%.3f",det[i,"z"] >posOut;
        }
    printf "\n" >posOut;
    
    for (i=1; i<=numDet; i++) {
        if ( ! det[i,"u"] ) continue;
        if (i!=1) printf "," >posOut;
        printf "%.3f",det[i,"th"]/DEG >posOut;
    }
    printf "\n" >posOut;
    
    for (i=1; i<=numDet; i++) {
        if ( ! det[i,"u"] ) continue;
        if (i!=1) printf "," >posOut;
        printf "%.3f",det[i,"ph"]/DEG >posOut;
    }
    printf "\n" >posOut;
}

function initGlobalVar (c) { # initialize global variables
# some constants
    PI = 4.*atan2(1.0,1.0);
    DEG = PI/180.;

# default beam velocity
    beta = 0.6;

# thickness of housing (added to each side of each crystal) in cm 
    wall = .1; 

# possible directions
    direct["x"] = 0;
    direct["y"] = 0;
    direct["z"] = 0;

# possible coordinates
    c["+"] = 0;
    c["0"] = 0;
    c["-"] = 0;
    for (x in c) {
        for (y in c) {
            for (z in c) {
                coord[x y z] = 0;  # corners and edge middle points
            }
        }
    }
    coord["c"]  = 0; # crystal center

# directions
    coord["xd"] = 0; # intrinsic x-direction of crystal
    coord["yd"] = 0; # intrinsic y-direction of crystal
    coord["zd"] = 0; # intrinsic z-direction of crystal
}

#
#
#  //////////////////////////
#
#

BEGIN {
    print "Running...";

    initGlobalVar();

    # read detector type information
    setupTypeSpec();
    
# just testing begin 
#    resetDet(1)
#    initDet(1,501)
#    placeDet(1,0,0,0,0)
#    transDet(1,"z",1)
#    rotDet(1,"x",45*DEG)
#    showDet(1);
#    exit
# just testing end
    

# start parsing the command line    
    cmdLine = "";
    for (i = 1; i < ARGC; i++) cmdLine = cmdLine " " ARGV[i];
    numStep = 0;
    for (i = 1; i < ARGC; i++) {
        
        if ( 0 ) {
        } else if ( ARGV[i] == "-a" ) {        # angle
            numStep++;
            ang[numStep,"l"] = ARGV[++i]*DEG;  # (lower) limit
            ang[numStep,"r"] = ARGV[++i]*DEG;  # resolution (command line given in degree)  or in %FWHM (for neg value)

        } else if ( ARGV[i] == "-b" ) {        # beta
            beta = ARGV[++i];                  

        } else if ( ARGV[i] == "-w" ) {        # wall thickness
            wall = ARGV[++i];                  

        } else if ( ARGV[i] == "-g" ) {        # group
            numGrpStep++;
            grp[numGrpStep,"l"] = ARGV[++i]*DEG;  # (lower) limit
            grp[numGrpStep,"n"] = ARGV[++i];   # integer number of detectors to group

        } else if ( ARGV[i] == "-t" ) {        # detector types
            numType++;
            type[numType,"t"] = ARGV[++i];     # detector type
            type[numType,"n"] = ARGV[++i];     # number of detectors of this type

        } else if ( ARGV[i] == "-v" ) deb = 1;
        else {
            print "Unknown option i, ARGV[i]:", i, ARGV[i];
            exit 1;
        }
    }

# finish processing of command line arguments
    ARGC=1;  # otherwise arguments are interpreted as input files later
    
    if ( ! numStep ) {  # default values for lower limit and resolution
        numStep++;
        ang[numStep,"l"] = 10*DEG;
        ang[numStep,"r"] = -10*DEG;
    }
    numStep++;
    ang[numStep,"l"] = 180*DEG;
    ang[numStep,"r"] = 0;

    if ( ! numGrpStep ) {  # default: do not group (grouping == 1)
        numGrpStep++;
        grp[numGrpStep,"l"] = 1*DEG;
        grp[numGrpStep,"n"] = 1;
    }
    numGrpStep++;
    grp[numGrpStep,"l"] = 180*DEG;
    grp[numGrpStep,"n"] = 1;


    if ( ! numType ) {  # default if no types are given
#
# there are 50 type 3 detectors   DALI1
#           84 type 1  DALI2
#           90 type 2  DALI2
#
        numType++;
        type[numType,"t"] = 3
        type[numType,"n"] = 6
        numType++;
        type[numType,"t"] = 1
        type[numType,"n"] = 82
        numType++;
        type[numType,"t"] = 2
        type[numType,"n"] = 96
        numType++;
        type[numType,"t"] = 3
        type[numType,"n"] = 62-type[1,"n"]
    }

    for (i=1; i<=numStep; i++) {
        print "ang", i, ang[i,"l"]/DEG, ang[i,"r"]/DEG;
    }

    for (i=1; i<=numGrpStep; i++) {
        print "grp", i, grp[i,"l"]/DEG, grp[i,"n"];
    }

    for (i = 1; i<=numType; i++) {
        typeTotal[type[i,"t"]] += type[i,"n"];
        printf "Using %2d detectors of type %d.\n", type[i,"n"], type[i,"t"];
    }
    for (t in typeTotal) {
        printf "Using a total of %3d detectors of type %d.\n", typeTotal[t], t;
    }

# setup detectors (i.e. fill the array det[])
    detSetup();
    print "Done with setup.  Now positioning the detectors.";

# start with positioning
# algorithm: - terminate when all detectors are used, or when full coverage reached
#    for placement: find distance such that delta_theta is just smaller than ang[is,"r"] assuming theta=ang[is,"t"]
#                   find angle such that theta_low is larger than ang[is,"l"]
#                   calculate phi coverage and increase that until one extra detector fits in (optimizing the phi coverage)
#                   if there are not detectors of the same type left, skip the rest and use the next type
#                   calculate theta_high and repeat procedure with
#                          ang[is,"l"] = theta_high or doing is++, if theta_high > ang[is+1,"l"]
#
#    for (i=1; i<=numDet; i++) resetPos(i);

    id = 1;  # detector number
    ig = 1;  # number of detectors to group
    is = 1;  # angle step number
    enRes =0; 
    ll = ang[1,"l"];   # Lower Limit for detector edge
    while (id <= numDet && is <= numStep) {   # while there are still detectors and rings to fill
        print "New ring ", iRing+1, ".  Working on detector", id
        if (deb) print "id, is, ig", id, is, ig
        tc = 180;
        tl = 0;
        maxLoop = 100000;
        rr = ang[is,"r"];    # aimed angulaR Resolution  (or energy resolution if negative)
        gr = grp[ig,"n"];    # number of detectors to group parallel 
        osx1 = det[id,"sx1"];  # save original x-sizes
        osx2 = det[id,"sx2"];
        osy1 = det[id,"sy1"];  # save original y-sizes
        osy2 = det[id,"sy2"];
        osz  = det[id,"sz"];   # save original z-sizes
        xshift = osx1 > osx2 ? osx1 : osx2;  # the shift of the position in crystal x-position for grouped crystals
                                             # should correspond to the largets size in x
        # create dummy detector (multiply x-size) if grouping is used
        det[id,"sx1"] *= gr;
        det[id,"sx2"] *= gr;

        # add wall thickness to each side
        det[id,"sx1"] += 2.0*wall
        det[id,"sx2"] += 2.0*wall
        det[id,"sy1"] += 2.0*wall
        det[id,"sy2"] += 2.0*wall
        det[id,"sz"] += 2.0*wall

        if (deb) print "start of new ring loop: id, rr, gr, osx1, osx2, sx1, sx2, xshift", \
                                                id, rr, gr, osx1, osx2, det[id,"sx1"], det[id,"sx2"], xshift;
        if (rr < 0) {   # ang[is,"r"] is the energy resolution in %FWHM, convert to angle
            ang[is,"l"]
            enRes = - 0.01 * rr/DEG;
            enRes *= gr;
            if (deb) print "enRes,ll/DEG", enRes, ll/DEG;
            rr = enRes * (1. - beta*cos(ll)) / (beta*sin(ll));  # start with lower angle and then iterate a couple of times
            if (deb) print "ll, rr", ll/DEG, rr/DEG;
            rr = enRes * (1. - beta*cos(ll+rr/2)) / (beta*sin(ll+rr/2));
            if (deb) print "ll, rr", ll/DEG, rr/DEG;
            rr = enRes * (1. - beta*cos(ll+rr/2)) / (beta*sin(ll+rr/2));
            if (deb) print "ll, rr", ll/DEG, rr/DEG;
            rr = enRes * (1. - beta*cos(ll+rr/2)) / (beta*sin(ll+rr/2));
            if (deb) print "ll, rr", ll/DEG, rr/DEG;
            rr = enRes * (1. - beta*cos(ll+rr/2)) / (beta*sin(ll+rr/2));
            if (deb) print "ll, rr", ll/DEG, rr/DEG;
            rr = enRes * (1. - beta*cos(ll+rr/2)) / (beta*sin(ll+rr/2));
            if (deb) print "ll, rr", ll/DEG, rr/DEG;
            rr = enRes * (1. - beta*cos(ll+rr/2)) / (beta*sin(ll+rr/2));
        }
        if (deb) print "ll, rr", ll/DEG, rr/DEG;
        rrOrig=rr;
        if (ll + rr > 175*DEG) rr = 175*DEG-ll;     # upper placement limit of 175 degree
        if (rr < rrOrig/2) break;                   # no space left for another ring

        print "Aiming for angular resolution of ", rr/DEG, " deg.";
        if (deb) print "rr", rr/DEG;
        d =  det[id,"sx2"]/2./tan(rr/2.);           # starting value for distance (based on front face of detector)
        d += det[id,"sz"]/2.;
        t = ll + rr/2.;                             # starting values for placement angle (Theta of detector center)
                                                    # lower limit + 1/2 of the aimed resolution
        if (deb) printf "    ";
        enReCalc = 2;  # number of recalculations of angular resolution (only needed if energy resolution is given)
                       # recalculation, since initial theta value might have changed considerably during placement
        if ( ! enRes ) enReCalc = 1;
        while (enReCalc--) {
            rrOld = rr;
            if ( enRes ) {
                rr    = enRes * (1. - beta*cos(t)) / (beta*sin(t));
                if (rrOld < rr) rr = rrOld;  # was already good enough
                if (deb) print "After energy recalc ", enReCalc, " ll, rr, rrOld ", ll/DEG, rr/DEG, rrOld/DEG
            }
            while ( 1 ) {   # placement loop:   the loop is broken if one ring is nicely filled
                iPlacementSteps++;
                while ( t < 0  ) t += PI;
                while ( t > PI ) t -= PI;
                if ( d < 0 ) { print "ERR: d<0!!!"; exit 1 }
                resetPos(id);
                placeDet(id, d, t);    # place detector and determine the maximum and minimum theta and the theta and phy coverage
                thetaPhiCover(id,"dummy");
                tl = det[id,"tp-tl"];   # theta low
                th = det[id,"tp-th"];   # theta high
                tc = det[id,"tp-tc"];   # theta cover
                pc = det[id,"tp-pc"];   # phi cover
                pr = 2*PI % pc;   # phi remainder
                if (deb) printf "d=%7.3f t=%8.3f tc=%6.2f (%7.3f, rr=%7.3f) tl=%6.2f (ll=%6.2f) pc=%6.2f pr=%6.2f 2pi/pc=%6.2f\n", \
                             d,      t/DEG,  tc/DEG, ang[is,"r"]/DEG, rr/DEG, tl/DEG, ll/DEG, pc/DEG, pr/DEG, 2.*PI/pc;
                if (! maxLoop--) { print "MAXLOOP!!!!"; exit 1; }
                if ( tl > ll + 0.1*DEG ) {   # lower theta of detector is too large
                    if (deb) printf "tl> ";
                    t -= tl - (ll + 0.03*DEG);
                    continue;
                }
                if ( tc < rr - 0.1*DEG ) { # theta limit overstepped by too much
                    if (deb) printf "tc< ";
                    del = (d - det[i,"sz"]/2.) * ( rr/tc - 1. );
                    d -= del > 0.001 ? del : 0.001;
                    continue;
                }
                if ( tc > rr           ) { # theta coverage too large; slowly increase d until limit reached
                    if (deb) printf "tc> ";
                    del = (d - det[i,"sz"]/2.) * ( tc/rr - 1. );
                    del *= 0.5*sin(t);
                    d += del > 0.001 ? del : 0.001;
                    if (del> 0.001) t -= 0.4*(tc-rr)/2.;
                    continue;
                }
                if ( tl < ll           ) { # lower theta is too small; must increase t until limit reached
                    if (deb) printf "tl< ";
                    t += ( ll - tl ) + 0.01*DEG;
                    continue;
                }
                if ( pr > 1*DEG        ) { # phi coverage not reached; slowry increase wanted resolution 
                                           # until another detector fits in
                    if (deb) printf "pr> ";
                    del = rr * 0.5 * (pc-pr)/(2*PI)*sin(t);
                    rr -= del > 0.00001 ? del : 0.00001;
                    continue;
                }
                break;
            }  # while (1) # placement loop
        }  # while enReCalc

        if (deb) printf "d, t, tc, tl, pc, pr %7.3f %6.1f %6.2f (%6.2f, %6.2f) %6.2f (%6.2f) %6.2f %6.2f\n", \
                     d, t/DEG, tc/DEG, ang[is,"r"]/DEG, rr/DEG, tl/DEG, ll/DEG, pc/DEG, pr/DEG;
        n = int(360*DEG/pc);

        printf "%s", "Need to place " n*gr " detectors.  ";
        if (id + n*gr - 1 > numDet) {
            print "Don't have enough detectors to fill ring", iRing, "!!!"
            n = numDet - id + 1;
        }
# check if there are n detectors of the same type available
        if ( det[id,"t"] != det[id + (n * gr) - 1,"t"] ) {
            print "Don't have enough detectors.";
            ii = id;
            while ( det[ii,"t"] == det[id,"t"] && id <= numDet ) {
                print "Skipping id", id, ".";
                detSkip[det[id,"t"]]++;
                id++;
            }
            continue;
        }

        # reset detector size
        det[id,"sx1"] = osx1;
        det[id,"sx2"] = osx2;
        det[id,"sy1"] = osy1;
        det[id,"sy2"] = osy2;
        det[id,"sz"]  = osz;

        ii = id + n*gr;
        for (i=1; i<=gr; i++){
            iRing++;
            x0 = -0.5*xshift*(gr-1) + xshift*(i-1);
            resetPos(id);
            placeDet(id,d,t,0,x0);
            thetaPhiCover(id,"dummy");

            numRing[iRing] = n;
            tRing[iRing]   = det[id,"000","theta"];
            tcRing[iRing]  = det[id,"tp-tc"];
            thRing[iRing]  = det[id,"tp-th"];
            tlRing[iRing]  = det[id,"tp-tl"];
            scRing[iRing]  = 0.5 * (cos(det[id,"tp-tl"]) - cos(det[id,"tp-th"]));  # solid angle coverage (= solidangle / 4pi)
            print "scring: ", scRing[iRing]

            print "Placing ", n, " detectors in ring", iRing, " theta coverage (",tl/DEG,",",th/DEG,"). Used ",iPlacementSteps," steps so far.";
            iii = id + n;
            while ( id < iii && id <= numDet ) {
                det[id,"r"] = iRing;
                det[id,"u"] = 1;
                det[id,"d"] = d;
                det[id,"x0"] = x0;
                det[id,"th"] = t;
                det[id,"ph"] = pc * ( id - (iii - n) );
                if (deb) print "final values: iRing, d, x0, th, ph", iRing, d,  det[id,"x0"], det[id,"th"]/DEG, det[id,"ph"]/DEG ;
                id++;
            }
        }

# update of lower limit angle
        if ( th > ang[is + 1,"l"] ) is++;
        ll = th;

# update grouping
        if ( th > grp[ig + 1,"l"] ) ig++;

#    if ( ang[is,"r"] + th > PI) {
        if ( rr + th > PI) {
            print "Array full.  Skipping detectors ", id, " to ", numDet;
            while (id <= numDet) detSkip[det[id++,"t"]]++;
            break;
        }

    }  #   while (id <= numDet && is <= numStep)   # while there are still detectors and ring to fill
    numOfRings = iRing;

#
#
# the main things are finished here
#  final placement and output from here on
#

    print "";
    print "# Output of angles (units: cm, degree)";
#  print "i | x y z | xd1-3 | yd1-3 | zd1-3 | theta phi";
    print "#          |           ----  center coordinates   ---- | -----   crystal size   -----";
    print "# det ring type  x     y      z  |   r    theta    phi |  sx1    sx2    sy1    sy2    sz";
    for (i=1; i<=numDet; i++) {
        if ( ! det[i,"u"] ) continue;
        resetPos(i);
        placeDet(i, det[i,"d"], det[i,"th"], det[i,"ph"], det[i,"x0"]);
        thetaPhiCover(i,"dummy");
        printf "%3d %3d %4d %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n", \
            i, det[i,"r"], det[i,"t"], det[i,"c","x"], det[i,"c","y"], det[i,"c","z"], \
            det[i,"d"], det[i,"th"]/DEG, det[i,"ph"]/DEG, \
            det[i,"sx1"], det[i,"sx2"], det[i,"sy1"], det[i,"sy2"], det[i,"sz"];
# calculate solid angle
        a = det[i,"sx1"] * det[i,"sy1"];
        A = 4*PI*(det[i,"d"] + det[i,"sz"]/2)**2;
        r = a/A;
        bf = sqrt(1.-beta**2)/(1. - beta*cos(det[i,"th"]));  # this is E_lab/E_cm
        bf = bf*bf;                                          # this is d_Omega_cm/d_Omega_lab
        if ( boostFactor[det[i,"r"]] == "" ) boostFactor[det[i,"r"]] = bf;
        solidAngle["t"]               += r;   # this is actually the coverage ratio
        solidAngle[det[i,"r"]]        += r;   # this is actually the coverage ratio
        solidAngleBoosted["t"]        += r*bf;
        solidAngleBoosted[det[i,"r"]] += r*bf;
#    print a, A, a/A, solidAngle["t"], bf, solidAngleBoosted["t"];
    }

    if (0) {  # output of ring coverage
        for (i = 1; i <= numOfRings; i++) {
            t = tRing[i];
            printf "ring %3d %7.3f %7.3f %7.3f %7.3f %7.3f    %6.2f %6.4f %6.2f\n", i, solidAngle[i], solidAngleBoosted[i], boostFactor[i], scRing[i], solidAngle[i]/scRing[i], t/DEG, beta*sin(t)/(1. - beta*cos(t))*tcRing[i], tcRing[i]/DEG;
        }
        print "ring coverage tot", solidAngle["t"], solidAngleBoosted["t"]
    }

#
# gnuplot theta
    gpf = "t.gp_";  # gnuplot file
    print "set term post enhanced color solid \"Times-Roman\" 14" > gpf;
    print "set size ratio -1" >gpf;
    print "set multiplot" >gpf;
    print "set size .59,.6" >gpf;
    print "set rmargin 0; set lmargin 0; set tmargin 0; set bmargin 0" >gpf;
    print "set origin .1,.1" >gpf;
#  print "set xrange [-50:100]" >gpf;
    print "set yrange [0:*]" >gpf;
    print "" >gpf;
    print "set label \"{/=24 DALI 3 Configuration \\\"Sea Urchin\\\"}\" at screen .5,.99 cent" >gpf;
    print "set label \"{/Courier=10 Command line: ", cmdLine, "}\" at screen .1,.95" > gpf;
    printf "set label \"{/=10 solid angle coverage (boosted " beta "): %.2f (%.2f)\" at graph .95,.95 right\n", \
        solidAngle["t"], solidAngleBoosted["t"] > gpf;
    print "set label \"{/=10 Numbers inside the detectors: ring#: {/Symbol q}-coverage, det-type}\" at screen .1,.03" > gpf;
    print "set label \"{/=10 Numbers inside the ring: number of detectors in ring\" at screen .1,.01" > gpf;
    print "set label \"{/=10 Red things: input parameters (from command line)\" at screen .1,-.01 tc lt 1" > gpf;
    xpos = .1;
    ypos = .93;
    for (t in detSkip) {
        print "set label \"{/tt=10 Skipped ", detSkip[t], "detectors of type ", t, "}\" at screen " xpos "," ypos > gpf;
        ypos -= .02;
    }
    for (i=1; i<=numStep-1; i++) {
        r = 17;
        t = ang[i,"l"];
        print "set arrow from 0,0 to ",r*cos(t), ",", r*sin(t), "lt 1 lw 2 nohead" > gpf;
        print "set label \"{/=10 " t/DEG, "," ang[i,"r"]/DEG "}\" at ",r*cos(t), ",", r*sin(t) ,"rot by ", t/DEG " tc lt 1" >gpf;
    }
    for (i=1; i<=numGrpStep-1; i++) {
        r = 14;
        t = grp[i,"l"];
        print "set arrow from 0,0 to ",r*cos(t), ",", r*sin(t), "lt 2 lw 2 nohead" > gpf;
        print "set label \"{/=10 " t/DEG, "," grp[i,"n"] "}\" at ",r*cos(t), ",", r*sin(t) ,"rot by ", t/DEG " tc lt 2" >gpf;
    }
    for (i=1; i<=numDet; i++) {
        if ( ! det[i,"u"] ) continue;
        if ( det[i,"ph"] != 0 ) continue;
        ro  = det[i,"d"] + det[i,"sz"]/2;
        r  = det[i,"d"];
        t  = det[i,"000","theta"];
        tl = det[i,"tp-tl"];
        tc = det[i,"tp-tc"];
        print  "set arrow from 0,0 to " ro*cos(tl) "," ro*sin(tl) " lt 7 lw .5 nohead" >gpf;
        print  "set label \"{/=10 " tl/DEG, "}\" at ",ro*cos(tl), ",", ro*sin(tl) ," rot by ", tl/DEG >gpf;
        printf "set label \"{/=10 %d: %.2f, %d}\"%s\n",det[i,"r"], tc/DEG, det[i,"t"], " at " r*cos(t) "," r*sin(t) " cent rot by " t/DEG >gpf;
    }
    print "plot '-' notit w l lt 7 lw 2" >gpf;
    for (i=1; i<=numDet; i++) {
        if ( ! det[i,"u"] ) continue;
        if ( det[i,"ph"] != 0 ) continue;
        print det[i,"tp-shape xz"] > gpf;
        print "" >gpf;
    }
    print "e" >gpf;
#
# gnuplot phi
    xpos = .7;
    ypos = .1;
    xsiz = .1;
    ysiz = .1;
    print "unset label; unset arrow" >gpf;
    print "set size ratio -1" >gpf;
    print "set yrange [*:*]" >gpf;
    print "set size ",xsiz*.95,",",ysiz*.95 >gpf;
    print "set border 0 lw .5" >gpf;
#  print "set arrow from graph 0,0 rto 10,0 lc 7 lw 2 nohead" >gpf;   # rto does not seem to work
    print "unset xtics; unset ytics" >gpf;
    id = 1;
    ii = id;
    iRing=1;
    while (id <= numDet) {
        if ( ! det[id,"u"] ) {
            id++;
            continue;
        }
        ii = id;
        print "set origin ",xpos,",",ypos >gpf;
        print "set label 11 \"" iRing  "\" at graph 0.07,.95 right" >gpf;
        print "set label 12 \""  numRing[iRing++] "\" at 0,0 cent" >gpf;
        print "plot '-' notit w l lt 7" >gpf;
        while ( det[ii,"r"] == det[id,"r"] ) {
            print det[id,"tp-shape xy"] > gpf;
            id++;
        }
        print "e" > gpf;
        ypos += ysiz;
        if (ypos + ysiz >.905) {
            ypos = .1;
            xpos += xsiz;
        }
    }

    print "set size .45,.15; set size noratio;" > gpf;
    print "set origin .1,.75" > gpf;
    print "set border 15" > gpf;
    print "set xtics; set ytics nomirror; set y2tics nomirror" > gpf;
    print "set xrange [0:*]; set yrange [0:*]; set y2range [0:*]" > gpf;
    print "unset label" > gpf;
    print "set format x \"{/=10 %g}\"" > gpf;
    print "set format y \"{/=10 %g}\"" > gpf;
    print "set format y2 \"{/=10 %g}\"" > gpf;
    print "set label \"{/=8 Solid ang. coverage}\" at graph -.12,.5 rot cent" > gpf;
    print "set label \"{/=8 Resolution (FWHM)}\" at graph 1.1,.5 rot cent point lt 7 pt 7 ps .5" > gpf;
    print "plot '-' notit w l lt 1, '-' notit w l lt 7, '-' notit w lp lt 2 pt 1, '-' notit w lp lt 3 pt 1, '-' notit axis x1y2 w p lt 7 pt 7 ps .5" > gpf;
    for (i = 1; i <= numOfRings; i++) {
        print tlRing[i]/DEG,solidAngle[i]/scRing[i] >gpf;
        print thRing[i]/DEG,solidAngle[i]/scRing[i] >gpf;
        print > gpf;
    }
    print tlRing[1]/DEG,         solidAngle["t"] >gpf;
    print thRing[numOfRings]/DEG,solidAngle["t"] >gpf;
    print "e" > gpf;
    print tlRing[1]/DEG,         solidAngleBoosted["t"] >gpf;
    print thRing[numOfRings]/DEG,solidAngleBoosted["t"] >gpf;
    print "e" > gpf;
    for (i = 1; i <= numOfRings; i++) {
        t = tRing[i];
        print t/DEG, scRing[i]*2 >gpf;
    }
    print "e" > gpf;
    for (i = 1; i <= numOfRings; i++) {
        t = tRing[i];
        print t/DEG, scRing[i]*boostFactor[i]*2 >gpf;
    }
    print "e" > gpf;
    for (i = 1; i <= numOfRings; i++) {
        t = tRing[i];
        print t/DEG, beta*sin(t)/(1. - beta*cos(t))*tcRing[i] >gpf;
    }
    print "e" > gpf;

    system("mv t.gp_ t.gp");
    print "Run the following to get eps file"
    print "# gnuplot t.gp >t.eps"

#    cvsOutPut();

}
