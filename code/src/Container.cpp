#include "Container.h"

BOX::BOX() {
    for (int i = 0; i < 3; i++)
    {
        xboundary[i] = '\0';
        yboundary[i] = '\0';
        zboundary[i] = '\0';
    }
    xlo = xhi = 0;
    ylo = yhi = 0;
    zlo = zhi = 0;
    xy = xz = yz = 0;
    xlo_bound = xhi_bound = 0;
    ylo_bound = yhi_bound = 0;
    zlo_bound = zhi_bound = 0;
    triclinic = false;
    lx = ly = lz = vol = 0;
    for (int i = 0; i < 6; i++) {
        h[i] = 0;
        h_inv[i] = 0;
    }
    d[0] = d[1] = d[2] = 0;
}

BOX::~BOX() {

}

void BOX::set_box(bool bound) {
    double aXb[3], bXc[3], aXc[3];  // used for triclinic box. The cross product of two of three edge vectors: a=(lx,0,0), b=(xy,ly,0), c=(xz,yz,lz)
    double bound_array_1[4] = {0, xy, xz, xy + xz};
    double bound_array_2[2] = {0, yz};
    if (!triclinic && (xy != 0 || xz != 0 || yz != 0)) {
        fprintf(stderr, "The tilt factors of non-triclinic box should be 0 Error!\n");
        exit(1);
    }

    if (bound) {
        xlo = xlo_bound - *std::min_element(bound_array_1, bound_array_1 + 4);
        xhi = xhi_bound - *std::max_element(bound_array_1, bound_array_1 + 4);
        ylo = ylo_bound - *std::min_element(bound_array_2, bound_array_2 + 2);
        yhi = yhi_bound - *std::max_element(bound_array_2, bound_array_2 + 2);
        zlo = zlo_bound;
        zhi = zhi_bound;
    } else {
        xlo_bound = xlo + *std::min_element(bound_array_1, bound_array_1 + 4);
        xhi_bound = xhi + *std::max_element(bound_array_1, bound_array_1 + 4);
        ylo_bound = ylo + *std::min_element(bound_array_2, bound_array_2 + 2);
        yhi_bound = yhi + *std::max_element(bound_array_2, bound_array_2 + 2);
        zlo_bound = zlo;
        zhi_bound = zhi;
    }
    lx = xhi - xlo;
    ly = yhi - ylo;
    lz = zhi - zlo;
    vol = lx * ly * lz;
    h[0] = lx;
    h[1] = ly;
    h[2] = lz;
    h_inv[0] = 1.0 / h[0];
    h_inv[1] = 1.0 / h[1];
    h_inv[2] = 1.0 / h[2];
    if (triclinic) {
        h[3] = yz;
        h[4] = xz;
        h[5] = xy;
        h_inv[3] = -h[3] / h[1] / h[2];
        h_inv[4] = (h[3] * h[5] - h[1] * h[4]) / (h[0] * h[1] * h[2]);
        h_inv[5] = -h[5] / h[0] / h[1];
        aXb[0] = aXb[1] = 0;
        aXb[2] = lx * ly;
        bXc[0] = ly * lz;
        bXc[1] = -xy * lz;
        bXc[2] = xy * yz - xz * ly;
        aXc[0] = 0;
        aXc[1] = -lx * lz;
        aXc[2] = lx * yz;
        d[0] = vol / std::sqrt(bXc[0]*bXc[0] + bXc[1]*bXc[1] + bXc[2]*bXc[2]);
        d[1] = vol / std::sqrt(aXc[0]*aXc[0] + aXc[1]*aXc[1] + aXc[2]*aXc[2]);
        d[2] = vol / std::sqrt(aXb[0]*aXb[0] + aXb[1]*aXb[1] + aXb[2]*aXb[2]);
    } else {
        h[3] = h[4] = h[5] = 0;
        h_inv[3] = h_inv[4] = h_inv[5] = 0;
    }
}

STEP::STEP()
{
    Timestep = -1;
    NumOfAtoms = 0;
    Box = BOX();
    vector<double>().swap(Masses);
    vector<string>().swap(Properties);
    vector< vector<double> >().swap(Atoms);

    wrapped_flag = -1;
    scaled_flag = -1;
    ID=TYPE=MOL=DIAMETER=DENSITY=X=Y=Z=IX=IY=IZ=VX=VY=VZ=WX=WY=WZ=FX=FY=FZ=Q=-1;

    NumOfAtomTypes = 0;
    NumOfBonds = NumOfBondTypes = 0;
    NumOfAngles = NumOfAngleTypes = 0;
    NumOfDihedrals = NumOfDihedralTypes = 0;
    NumOfImpropers = NumOfImproperTypes = 0;
    vector< vector<int> >().swap(Bonds);
    vector< vector<int> >().swap(Angles);
    vector< vector<int> >().swap(Dihedrals);
    vector< vector<int> >().swap(Impropers);
}

STEP::~STEP()
{
    vector<double>().swap(Masses);
    vector<string>().swap(Properties);
    vector< vector<double> >().swap(Atoms);

    vector< vector<int> >().swap(Bonds);
    vector< vector<int> >().swap(Angles);
    vector< vector<int> >().swap(Dihedrals);
    vector< vector<int> >().swap(Impropers);
}

void STEP::x2lamda() {
    // convert Cartesian Coordinate to Scaled (lamda) Coordinate (i.e., using the three edge vectors of the simulation box as the basis)
    if ( X == -1 || Y == -1 || Z == -1 ) {
        fprintf(stderr, "Cannot find complete (X,Y,Z) properties!\n");
        return;
    }
    
    double dx, dy, dz;

    for (int i = 0; i < (int)Atoms.size(); i++) {
        dx = Atoms[i][X] - Box.xlo;
        dy = Atoms[i][Y] - Box.ylo;
        dz = Atoms[i][Z] - Box.zlo;

        Atoms[i][X] = Box.h_inv[0] * dx + Box.h_inv[5] * dy + Box.h_inv[4] * dz;
        Atoms[i][Y] = Box.h_inv[1] * dy + Box.h_inv[3] * dz;
        Atoms[i][Z] = Box.h_inv[2] * dz;
    }
}

void STEP::lamda2x() {
    // For triclinic box, convert lamda coordinate (coordinate udner parallelepiped base vectors) to Cartesian coordinate
    // x = H lamda + x0
    if ( X == -1 || Y == -1 || Z == -1 ) {
        fprintf(stderr, "Cannot find complete (X,Y,Z) properties!\n");
        return;
    }

    for (int i = 0; i < (int)Atoms.size(); i++) {
        Atoms[i][X] = Box.h[0] * Atoms[i][X] + Box.h[5] * Atoms[i][Y] + Box.h[4] * Atoms[i][Z] + Box.xlo;
        Atoms[i][Y] = Box.h[1] * Atoms[i][Y] + Box.h[3] * Atoms[i][Z] + Box.ylo;
        Atoms[i][Z] = Box.h[2] * Atoms[i][Z] + Box.zlo;
    }
}

int STEP::wrap() {
    bool no_image, with_image;

    if ( X == -1 || Y == -1 || Z == -1 ) {
        fprintf(stderr, "Cannot find complete (X,Y,Z) properties!\n");
        return 1;
    }
    no_image = (IX == -1 && IY == -1 && IZ == -1);
    with_image = (IX != -1 && IY != -1 && IZ != -1);

    if (scaled_flag == 0 && no_image) {      // unscaled (wrapped or unwrapped), without image flag
        x2lamda();
        for (int i = 0; i < (int)Atoms.size(); i++) {
            Atoms[i][X] -= floor(Atoms[i][X]);
            Atoms[i][Y] -= floor(Atoms[i][Y]);
            Atoms[i][Z] -= floor(Atoms[i][Z]);
        }
        lamda2x();
    } else if (scaled_flag == 0 && wrapped_flag == 0 && with_image) {  // unscaled and unwrapped, with image flag
        x2lamda();
        for (int i = 0; i < (int)Atoms.size(); i++) {
            Atoms[i][IX] = floor(Atoms[i][X]);
            Atoms[i][IY] = floor(Atoms[i][Y]);
            Atoms[i][IZ] = floor(Atoms[i][Z]);
            Atoms[i][X] -= floor(Atoms[i][X]);
            Atoms[i][Y] -= floor(Atoms[i][Y]);
            Atoms[i][Z] -= floor(Atoms[i][Z]);
        }
        lamda2x();
    } else if (scaled_flag == 0 && wrapped_flag > 0 && with_image) { // unscaled and wrapped, with image flag
        x2lamda();
        for (int i = 0; i < (int)Atoms.size(); i++) {
            Atoms[i][IX] += floor(Atoms[i][X]);
            Atoms[i][IY] += floor(Atoms[i][Y]);
            Atoms[i][IZ] += floor(Atoms[i][Z]);
            Atoms[i][X] -= floor(Atoms[i][X]);
            Atoms[i][Y] -= floor(Atoms[i][Y]);
            Atoms[i][Z] -= floor(Atoms[i][Z]);
        }
        lamda2x();
    } else if (scaled_flag == 1 && no_image) { // scaled, without image flag
        for (int i = 0; i < (int)Atoms.size(); i++) {
            Atoms[i][X] -= floor(Atoms[i][X]);
            Atoms[i][Y] -= floor(Atoms[i][Y]);
            Atoms[i][Z] -= floor(Atoms[i][Z]);
        }
        lamda2x();
    } else if (scaled_flag == 1 && wrapped_flag == 0 && with_image) { // scaled and unwrapped, with image flag
        for (int i = 0; i < (int)Atoms.size(); i++) {
            Atoms[i][IX] = floor(Atoms[i][X]);
            Atoms[i][IY] = floor(Atoms[i][Y]);
            Atoms[i][IZ] = floor(Atoms[i][Z]);
            Atoms[i][X] -= floor(Atoms[i][X]);
            Atoms[i][Y] -= floor(Atoms[i][Y]);
            Atoms[i][Z] -= floor(Atoms[i][Z]);
        }
        lamda2x();
    } else if (scaled_flag == 1 && wrapped_flag > 0 && with_image) { // scaled and wrapped, with image flag
        for (int i = 0; i < (int)Atoms.size(); i++) {
            Atoms[i][IX] += floor(Atoms[i][X]);
            Atoms[i][IY] += floor(Atoms[i][Y]);
            Atoms[i][IZ] += floor(Atoms[i][Z]);
            Atoms[i][X] -= floor(Atoms[i][X]);
            Atoms[i][Y] -= floor(Atoms[i][Y]);
            Atoms[i][Z] -= floor(Atoms[i][Z]);
        }
        lamda2x();
    }

    wrapped_flag = 1;
    scaled_flag = 0;
    for (int i = 0; i < (int)Properties.size(); i++)
        if (Properties[i] == "xu" || Properties[i] == "xs" || Properties[i] == "xsu" )
            Properties[i] = "x";
        else if (Properties[i] == "yu" || Properties[i] == "ys" || Properties[i] == "ysu")
            Properties[i] = "y";
        else if (Properties[i] == "zu" || Properties[i] == "zs" || Properties[i] == "zsu")
            Properties[i] = "z";

    return 0;
}

int STEP::unwrap() {
    if (X == -1 || Y == -1 || Z == -1) {
        fprintf(stderr, "Cannot find complete (X,Y,Z) properties!\n");
        return 1;
    }

    if (scaled_flag == 1 && wrapped_flag == 0) { // scaled and unwrapped
        lamda2x();
    } else if (scaled_flag == 0 && wrapped_flag == 0) {
    } else if (wrapped_flag > 0 ) { // wrapped_flag can be 1 or 2, both are wrapped (refer to explanation in header file)
        if (IX == -1 || IY == -1 || IZ == -1) {
            fprintf(stderr, "Cannot find complete (IX,IY,IZ) properties!\n");
            return 1;
        }

        if (scaled_flag == 1) { // scaled and wrapped
            for (int i = 0; i < (int)Atoms.size(); i++) {
                Atoms[i][X] += Atoms[i][IX];
                Atoms[i][Y] += Atoms[i][IY];
                Atoms[i][Z] += Atoms[i][IZ];
            }
            lamda2x();
        } else if (scaled_flag == 0) { // unscaled and wrapped
            x2lamda();
            for (int i = 0; i < (int)Atoms.size(); i++) {
                Atoms[i][X] += Atoms[i][IX];
                Atoms[i][Y] += Atoms[i][IY];
                Atoms[i][Z] += Atoms[i][IZ];
            }
            lamda2x();
        } else {
            fprintf(stderr, "Wrong Coordinate Type! Please check your input file's format!\n");
            return 1;
        }
    } else {
        fprintf(stderr, "Wrong Coordinate Type! Please check your input file's format!\n");
        return 1;
    }

    wrapped_flag = 0;
    scaled_flag = 0;
    for (int i = 0; i < (int)Properties.size(); i++)
        if (Properties[i] == "x" || Properties[i] == "xs" || Properties[i] == "xsu" )
            Properties[i] = "xu";
        else if (Properties[i] == "y" || Properties[i] == "ys" || Properties[i] == "ysu")
            Properties[i] = "yu";
        else if (Properties[i] == "z" || Properties[i] == "zs" || Properties[i] == "zsu")
            Properties[i] = "zu";

    return 0;
}

void STEP::flip() {
    // flip if tilt ratio > 0.5
    // NOTE: It assume that the tilt ratio is between -1 and 1.
    
    int xyflip = 0, xzflip = 0, yzflip = 0;
    // adust yz, xz and xy, the sequence cannot be changed
    if (Box.yz / Box.ly > 0.5) {
        Box.yz -= Box.ly;
        Box.xz -= Box.xy;
        yzflip--;
    } else if (Box.yz / Box.ly < -0.5) {
        Box.yz += Box.ly;
        Box.xz += Box.xy;
        yzflip++;
    }
    if (Box.xz / Box.lx > 0.5) {
        Box.xz -= Box.lx;
        xzflip--;
    } else if (Box.xz / Box.lx < -0.5) {
        Box.xz += Box.lx;
        xzflip++;
    }
    if (Box.xy / Box.lx > 0.5) {
        Box.xy -= Box.lx;
        xyflip--;
    } else if (Box.xy / Box.lx < -0.5) {
        Box.xy += Box.lx;
        xyflip++;
    }
    Box.set_box(false);

    /* ---------This is from LAMMPS 7Aug19-------------------------------------
    adjust image flags due to triclinic box flip (m=xyflip, n=xzflip, p=yzflip)
    flip operation is changing box vectors A,B,C to new A',B',C'
        A' = A              (A does not change)
        B' = B + mA         (B shifted by A)
        C' = C + pB + nA    (C shifted by B and/or A)
    this requires the image flags change from (a,b,c) to (a',b',c')
    so that x_unwrap for each atom is same before/after
        x_unwrap_before = xlocal + aA + bB + cC
        x_unwrap_after = xlocal + a'A' + b'B' + c'C'
    this requires:
        c' = c
        b' = b - cp
        a' = a - (b-cp)m - cn = a - b'm - cn
    in other words, for xy flip, change in x flag depends on current y flag
    this is b/c the xy flip dramatically changes which tiled image of
        simulation box an unwrapped point maps to
    ------------------------------------------------------------------------- */
    if (xyflip != 0 || xzflip != 0 || yzflip != 0) {
        #pragma omp parallel for
        for (int i = 0; i < NumOfAtoms; i++) {
            Atoms[i][IY] -= yzflip * Atoms[i][IZ];
            Atoms[i][IX] -= xyflip * Atoms[i][IY] + xzflip * Atoms[i][IZ];
        }
        this->wrap();
    }
}