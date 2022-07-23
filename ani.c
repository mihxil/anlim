/* ani.f -- translated by f2c (version 19970805).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Common Block Declarations */

struct {
    doublereal t2p;
    integer maxr, maxi;
} lookuc_;

#define lookuc_1 lookuc_

struct {
    doublereal rmul, radd, rmul2, radd2, rmul3, radd3;
} lookup_;

#define lookup_1 lookup_

union {
    struct {
	integer nrannr, irn[480000];
    } _1;
    struct {
	integer n, irn[480000];
    } _2;
} nran_;

#define nran_1 (nran_._1)
#define nran_2 (nran_._2)

struct {
    doublereal zsbnd[240000]	/* was [1600][150] */;
    shortint iz0s[1600], nin[1600];
} mats_;

#define mats_1 mats_

struct {
    shortint ngc[6400]	/* was [4][1600] */;
} ngcb_;

#define ngcb_1 ngcb_

struct {
    doublereal a3, wn;
    integer n1, n2, n12;
} dats_;

#define dats_1 dats_

struct {
    doublereal p[9], q[81]	/* was [9][9] */, pp[18000]	/* was [9][
	    2000] */;
    integer ipr[2000];
} pqrh_;

#define pqrh_1 pqrh_

struct {
    doublereal b[81]	/* was [9][9] */;
} cbh_;

#define cbh_1 cbh_

struct {
    integer ir1[9689], ir2[127];
    shortint inxt1[9689], inxt2[127];
} ransrb_;

#define ransrb_1 ransrb_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__5 = 5;
static real c_b37 = (float)0.;
static integer c_b98 = 480000;
static integer c__150 = 150;

/* Version: 1999-2-4 */

/*   MM: */
/*   hardly changed anything until now. Mainly added commence and lay out. */



/* for wolff cluster 3d s=1/2 ising model anisotropic limit */

/* select neighbors via continuous intervals to be skipped */
/* nrannr,irn in common block */
/* periodic boundary conditions */
/* program anlim.f (based on w1n.f) */

/* core: wn in common; vz=zsbnd ... */
/* a little bit cryptic (MM) */

/* Main program */ MAIN__()
{
    /* System generated locals */
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle(), s_rsle(), e_rsle();
    /* Subroutine */ int s_stop();
    integer f_open(), f_clos();

    /* Local variables */
    extern /* Subroutine */ int bain_();
    static doublereal tpar;
    static integer nint, itmax;
    static doublereal a3;
    static integer n1, n2, ntoss, nr, ncycle, iprlll, numint;
    static doublereal fac;

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, 0, 0 };
    static cilist io___2 = { 0, 5, 0, 0, 0 };
    static cilist io___12 = { 0, 6, 0, 0, 0 };
    static cilist io___13 = { 0, 5, 0, 0, 0 };


/* MM */
L102:
/*     read from standard input (MM) */
/* begin of eternal loop (MM) */
    s_wsle(&io___1);
    do_lio(&c__9, &c__1, " n1,n2,a3,ncycle,ntoss,numint,nint,nr,itmax", 43L);
    e_wsle();
    s_rsle(&io___2);
    do_lio(&c__3, &c__1, (char *)&n1, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&n2, (ftnlen)sizeof(integer));
    do_lio(&c__5, &c__1, (char *)&a3, (ftnlen)sizeof(doublereal));
    do_lio(&c__3, &c__1, (char *)&ncycle, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&ntoss, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&numint, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&nint, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&nr, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&itmax, (ftnlen)sizeof(integer));
    e_rsle();
    itmax *= 3600;
    if (n1 == 0) {
	s_stop("", 0L);
    }
/*     read more from standard input (MM) */
/* that's how we break out (MM) */
    s_wsle(&io___12);
    do_lio(&c__9, &c__1, " tpar,fac", 9L);
    e_wsle();
    s_rsle(&io___13);
    do_lio(&c__5, &c__1, (char *)&tpar, (ftnlen)sizeof(doublereal));
    do_lio(&c__5, &c__1, (char *)&fac, (ftnlen)sizeof(doublereal));
    e_rsle();
/*     open the output file (MM) */
    o__1.oerr = 0;
    o__1.ounit = 8;
    o__1.ofnmlen = 8;
    o__1.ofnm = "anli.dat";
    o__1.orl = 0;
    o__1.osta = 0;
    o__1.oacc = "append";
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
/*     main loop, which is not a loop right now, because nprlll = 1 (MM) 
*/
    for (iprlll = 1; iprlll <= 1; ++iprlll) {
	bain_(&n1, &n2, &a3, &ncycle, &ntoss, &numint, &nint, &nr, &itmax, &
		tpar, &fac);
/* L20: */
    }
/*     close the output file (MM) */
    cl__1.cerr = 0;
    cl__1.cunit = 8;
    cl__1.csta = 0;
    f_clos(&cl__1);
    goto L102;
/* eternal loop (MM) */
} /* MAIN__ */

/*     end program anlim (MM) */
/* ---------------------------------------------------------- */
/*  main program */
/*  (why it's bain, not main, I don't know.(MM) */
/* ---------------------------------------------------------- */
/* Subroutine */ int bain_(n1, n2, a3, ncycle, ntoss, numint, nint, nr, itmax,
	 tpar, fac)
integer *n1, *n2;
doublereal *a3;
integer *ncycle, *ntoss, *numint, *nint, *nr, *itmax;
doublereal *tpar, *fac;
{
    /* Initialized data */

    static char ident[80+1] = "anli v0                                      \
                                   ";

    /* Format strings */
    static char fmt_900[] = "(\002program for mc simulation  on a\002,i4,\
\002*\002,i4,\002*\002,f12.6,\002 lattice\002)";
    static char fmt_901[] = "(\002 tpar= \002,f11.6,\002 factor=\002,f11.6)";
    static char fmt_902[] = "(1x,i3,\002 consecutive runs of \002,i10,\002 W\
olff steps\002,1x,\002after\002,i10,\002 times nint steps for equilibratio\
n\002/)";
    static char fmt_924[] = "(\002 nr=\002,i8/)";
    static char fmt_903[] = "(\002 data are taken at intervals of \002,i3\
,\002 mcs/site\002/)";
    static char fmt_120[] = "(a8,\002; n1,numint,nint,knn,s:     \002,i3,i10\
,i3,f12.6,12x,i4)";
    static char fmt_554[] = "(20x,\002setup time (sec.)=\002,f10.6)";
    static char fmt_999[] = "(/70(\002*\002)/)";
    static char fmt_1120[] = "(\002 ****** emergency stop,mcs=\002,i9,\002\
 = \002,f11.3,\002 sec.\002)";
    static char fmt_1130[] = "(\002 mc stop at t=\002,f12.2,\002 sec.\002)";
    static char fmt_553[] = "(\002max nr random numbers used:\002,i10,\002 a\
llowed:\002,i10/\002max nr  cluster boundaries:\002,i10,\002 allowed:\002,i1\
0)";
    static char fmt_555[] = "(\0020mc simulation  time (microseconds/mc step\
/site)=\002,f20.6//\002 summation over lattice (microseconds/time/site)=\002\
,f20.6//\002 analysis (seconds) =\002,f10.2//\002 total time (sec.)=\002,f10\
.2)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsle(), do_lio(), e_wsle();
    /* Subroutine */ int s_stop();
    integer s_wsfe(), do_fio(), e_wsfe();

    /* Local variables */
    static doublereal anal;
    static integer ncup;
    static doublereal tmax;
    extern /* Subroutine */ int part_();
    extern doublereal ctijd_();
    extern /* Subroutine */ int carlo_();
    static integer nsigh, ntemp, npart;
    extern /* Subroutine */ int monte_();
    static integer itoss;
    static doublereal t4;
    static integer nc, ip;
    extern /* Subroutine */ int addpqr_();
    static integer ncount, mcstps;
    extern /* Subroutine */ int setpqr_();
    static doublereal xt0, xt1, xt2, xt3, xt4;
    static integer iip;
    static doublereal tan__, tmc;
    static integer nip, inp, npp;
    static doublereal vcu;
    extern /* Subroutine */ int pqr_();

    /* Fortran I/O blocks */
    static cilist io___18 = { 0, 6, 0, 0, 0 };
    static cilist io___22 = { 0, 6, 0, 0, 0 };
    static cilist io___29 = { 0, 6, 0, fmt_900, 0 };
    static cilist io___30 = { 0, 6, 0, fmt_901, 0 };
    static cilist io___31 = { 0, 6, 0, fmt_902, 0 };
    static cilist io___32 = { 0, 6, 0, fmt_924, 0 };
    static cilist io___33 = { 0, 6, 0, fmt_903, 0 };
    static cilist io___36 = { 0, 8, 0, fmt_120, 0 };
    static cilist io___39 = { 0, 6, 0, fmt_554, 0 };
    static cilist io___42 = { 0, 6, 0, fmt_999, 0 };
    static cilist io___49 = { 0, 6, 0, 0, 0 };
    static cilist io___50 = { 0, 6, 0, fmt_1120, 0 };
    static cilist io___51 = { 0, 6, 0, fmt_1130, 0 };
    static cilist io___52 = { 0, 6, 0, fmt_999, 0 };
    static cilist io___53 = { 0, 6, 0, fmt_553, 0 };
    static cilist io___57 = { 0, 6, 0, fmt_555, 0 };


/*     program for simulation of anisotropic limit of ising model */
    s_wsle(&io___18);
    do_lio(&c__9, &c__1, ident, 80L);
    e_wsle();
    xt0 = ctijd_(&c_b37);
    tmax = (doublereal) (*itmax - 10);
    ntemp = 0;
    if (*n1 > 40 || *n2 > 40) {
	s_wsle(&io___22);
	do_lio(&c__9, &c__1, " main: at least one lattice dimension too big", 
		45L);
	e_wsle();
	s_stop("", 0L);
    }
    mcstps = *numint * *nint;
    npart = *numint / 2;
    if (npart > 2000) {
	npart = 2000;
    }
    if (npart <= 10) {
	npart = 1;
    }
    npp = mcstps / npart;
    nip = npp / *nint;
    npp = *nint * nip;
    mcstps = npart * npp;
    ncount = mcstps / *nint;
/*     104   continue */
    ++ntemp;
/*     initialise: (MM) */
    carlo_(tpar, fac, n1, n2, a3, &ntemp, nr);
    vcu = *n1 * *n2 * *a3;
/*     write to standard output what we're doing (MM) */
    s_wsfe(&io___29);
    do_fio(&c__1, (char *)&(*n1), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*n2), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*a3), (ftnlen)sizeof(doublereal));
    e_wsfe();
/* added quotes and took */
    s_wsfe(&io___30);
    do_fio(&c__1, (char *)&(*tpar), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&(*fac), (ftnlen)sizeof(doublereal));
    e_wsfe();
    s_wsfe(&io___31);
    do_fio(&c__1, (char *)&(*ncycle), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&mcstps, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*ntoss), (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___32);
    do_fio(&c__1, (char *)&(*nr), (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___33);
    do_fio(&c__1, (char *)&(*nint), (ftnlen)sizeof(integer));
    e_wsfe();
    i__1 = *ntoss;
    for (itoss = 1; itoss <= i__1; ++itoss) {
	monte_(nint);
/* L6: */
    }
    i__1 = *ncycle;
    for (nsigh = 1; nsigh <= i__1; ++nsigh) {
	s_wsfe(&io___36);
	do_fio(&c__1, ident, 80L);
	do_fio(&c__1, (char *)&(*n1), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*numint), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*nint), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*tpar), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*nr), (ftnlen)sizeof(integer));
	e_wsfe();
	setpqr_();
	xt3 = ctijd_(&c_b37);
	xt1 = xt3 - xt0;
	s_wsfe(&io___39);
	do_fio(&c__1, (char *)&xt1, (ftnlen)sizeof(doublereal));
	e_wsfe();
	tmc = 0.;
	tan__ = 0.;
	s_wsfe(&io___42);
	e_wsfe();
	nc = 0;
	i__2 = npart;
	for (ip = 1; ip <= i__2; ++ip) {
	    inp = ip;
	    i__3 = nip;
	    for (iip = 1; iip <= i__3; ++iip) {
		++nc;
		monte_(nint);
		xt2 = ctijd_(&c_b37);
		tmc = tmc + xt2 - xt3;
		addpqr_(n1, n2, &ncup);
		xt3 = ctijd_(&c_b37);
		tan__ = tan__ + xt3 - xt2;
		if (tmax - xt3 + xt0 < 0.) {
		    goto L1111;
		}
/* L1: */
	    }
	    part_(&inp, &nc, &ncup);
/* L11: */
	}
	goto L1112;
L1111:
	ncount = nc;
	part_(&inp, &nc, &ncup);
	mcstps = ncount * *nint;
	s_wsle(&io___49);
	do_lio(&c__9, &c__1, "dat gaat niet goed", 18L);
	e_wsle();
	s_wsfe(&io___50);
	do_fio(&c__1, (char *)&mcstps, (ftnlen)sizeof(integer));
	d__1 = xt3 - xt0;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsfe();
L1112:
	s_wsfe(&io___51);
	d__1 = xt3 - xt0;
	do_fio(&c__1, (char *)&d__1, (ftnlen)sizeof(doublereal));
	e_wsfe();
/*     call lpsb */
	tmc = tmc * 1e6 / mcstps / vcu;
	tan__ = tan__ * 1e6 / ncount / vcu;
	pqr_(&ncup, &ncount, &npart, n1, n2, a3);
	if (tmax - xt3 + xt0 < 0.) {
	    s_stop("", 0L);
	}
/* L7: */
    }
    s_wsfe(&io___52);
    e_wsfe();
/* format(/70('*')/)   (MM) */
    s_wsfe(&io___53);
    do_fio(&c__1, (char *)&lookuc_1.maxr, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&c_b98, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&lookuc_1.maxi, (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&c__150, (ftnlen)sizeof(integer));
    e_wsfe();
    xt4 = ctijd_(&c_b37);
    anal = xt4 - xt3;
    t4 = xt4 - xt0;
    s_wsfe(&io___57);
    do_fio(&c__1, (char *)&tmc, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&tan__, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&anal, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&t4, (ftnlen)sizeof(doublereal));
    e_wsfe();
    return 0;
} /* bain_ */

/*     end bain (MM) */
/* ------------------------------------------------------ */
/* carlo */
/* Initialisation of Monte-Carlo algorithm.(MM) */

/* Subroutine */ int carlo_(tpar, fac, m1, m2, aa, ntmp, nr)
doublereal *tpar, *fac;
integer *m1, *m2;
doublereal *aa;
integer *ntmp, *nr;
{
    /* Format strings */
    static char fmt_9[] = "(\0021carlo disagrees with specified lattice siz\
e\002,2i3)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();
    /* Subroutine */ int s_stop();
    integer s_wsle(), do_lio(), e_wsle();

    /* Local variables */
    static integer itel, iseed, index;
    extern /* Subroutine */ int ransi_();
    static integer istep;
    extern /* Subroutine */ int ransr_();
    static integer ij, ix, iy;
    static doublereal rnd;
    static integer ipm, ixm, iym, ixp, iyp;
    static doublereal facc;

    /* Fortran I/O blocks */
    static cilist io___58 = { 0, 6, 0, fmt_9, 0 };
    static cilist io___69 = { 0, 6, 0, 0, 0 };


/*     version for periodic ising model */
/*     maximum sizes(MM): */
/*     check lattice size(MM): */
/*     in x-direction(MM) */
/* L4: */
    if (*m1 - 40 <= 0) {
	goto L8;
    } else {
	goto L7;
    }
/*     if wrong (for both x and y directions)(MM): */
L7:
    s_wsfe(&io___58);
    do_fio(&c__1, (char *)&(*m1), (ftnlen)sizeof(integer));
    do_fio(&c__1, (char *)&(*m2), (ftnlen)sizeof(integer));
    e_wsfe();
    s_stop("", 0L);
L8:
/*     in y-direction(MM) */
    if (*m2 - 40 <= 0) {
	goto L6;
    } else {
	goto L7;
    }
L6:
    dats_1.a3 = *aa;
    lookuc_1.t2p = *tpar / 2.;
    lookuc_1.maxr = 0;
    lookuc_1.maxi = 0;
    facc = *fac;
/*     size of system(MM) */
    dats_1.n1 = *m1;
    dats_1.n2 = *m2;
    dats_1.n12 = dats_1.n1 * dats_1.n2;
    dats_1.wn = 1. / (dats_1.a3 * dats_1.n12);
/*     istep is not used currently, I think (MM) */
    istep = 2;
/*     itel is probably some counter(MM) */
    itel = 0;
/*     This loop fills ngc. */
/*     this is a tabular in which the neigbors of every site are stored. 
*/
/*     (for the time being this is my theory for this (MM)) */
    i__1 = dats_1.n12;
    for (index = 1; index <= i__1; ++index) {
/*     write(6,202) index,itel,istep */
/*     02   format(' index,itel,istep',30i3) */
	++itel;
	iy = (index - 1) / dats_1.n1 + 1;
	ix = index - (iy - 1) * dats_1.n1;
	ixp = ix + 1 - ix / dats_1.n1 * dats_1.n1;
	iyp = iy + 1 - iy / dats_1.n2 * dats_1.n2;
	ixm = ix - 1 + (dats_1.n1 - ix + 1) / dats_1.n1 * dats_1.n1;
	iym = iy - 1 + (dats_1.n2 - iy + 1) / dats_1.n2 * dats_1.n2;
	ngcb_1.ngc[(index << 2) - 4] = (shortint) ((iy - 1) * dats_1.n1 + ixm)
		;
	ngcb_1.ngc[(index << 2) - 2] = (shortint) ((iy - 1) * dats_1.n1 + ixp)
		;
	ngcb_1.ngc[(index << 2) - 3] = (shortint) ((iym - 1) * dats_1.n1 + ix)
		;
	ngcb_1.ngc[(index << 2) - 1] = (shortint) ((iyp - 1) * dats_1.n1 + ix)
		;
/*     write(6,950)index,ngc(1,index),ngc(2,index),ngc(3,index), */
/*    ,ngc(4,index) */
/* 950  format('carlo index,ngc',i4,2x,4i4) */
/* L210: */
    }
/*     these things relate to the random-number generator:(MM) */
    lookup_1.radd = dats_1.n12 * .5 + 1.;
    lookup_1.rmul = dats_1.n12 * 2.328306436538694e-10;
    lookup_1.radd2 = dats_1.a3 * .5;
    lookup_1.rmul2 = dats_1.a3 * 2.328306436538694e-10;
    lookup_1.radd3 = .5;
    lookup_1.rmul3 = 2.328306436538694e-10;
    s_wsle(&io___69);
    do_lio(&c__9, &c__1, " wolff mc of s=1/2 anisotropic ising model", 42L);
    e_wsle();
    if (*ntmp > 1) {
	goto L20;
    }
/*     20 = jump to end of procedure. */
/*     'ntmp' probable means 'tempory integer' */
/*     it seems to be always .gt.1, exept for the first time.(MM) */
    iseed = abs(*nr) + 1;
/*     so 'nr' is what is called the seed from outside (MM) */
/*     initialize random number generator: */
    ransi_(&iseed);
    if (*nr >= 0) {
	goto L320;
    }
/*     if given 'seed' < 0 then we do this (MM): */
    i__1 = dats_1.n12;
    for (ij = 1; ij <= i__1; ++ij) {
	mats_1.iz0s[ij - 1] = 1;
	mats_1.nin[ij - 1] = 0;
/* L112: */
    }
/*     and we jump to end of this procedure(MM) */
/* MM */
    goto L20;
/*     nr >= 0: (MM) */
L320:
    nran_1.nrannr = dats_1.n12;
    ransr_();
    i__1 = dats_1.n12;
    for (ij = 1; ij <= i__1; ++ij) {
	rnd = nran_1.irn[ij - 1] * lookup_1.rmul3 + lookup_1.radd3;
	ipm = (integer) (rnd + rnd);
	ipm = ipm - 1 + ipm;
	mats_1.iz0s[ij - 1] = (shortint) ipm;
	mats_1.nin[ij - 1] = 0;
/* L330: */
    }
L20:
    nran_1.nrannr = dats_1.n12 * 300;
    return 0;
} /* carlo_ */

/*     end carlo (MM) */
/* ----------------------------------------------------- */
/* monte */
/* This subroutine only calls mcw (MM). */

/* Subroutine */ int monte_(nsteps)
integer *nsteps;
{
    extern /* Subroutine */ int mcw_();

/*     perform nsteps mc steps per site */
/*     version for arbitrary model */
/*     everything is implicitly a double exept i,j,k..n. */
/*     again the maximal sizes (MM): */
/*     do 2 nn=1,nsteps */
    mcw_(nsteps);
/*     2    continue */
    return 0;
} /* monte_ */

/* ---------------------------------------------------------- */
/* mcw  = monte */
/* This is the heart of the program, here happen */
/* the real Monte-Carlo simulations (MM) */
/* Subroutine */ int mcw_(nsteps)
integer *nsteps;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    double log();
    integer s_wsle(), do_lio(), e_wsle();
    /* Subroutine */ int s_stop();

    /* Local variables */
    static integer isnb, icsp;
    extern /* Subroutine */ int ransr_();
    static doublereal cl, dl;
    static integer ii, ij, ik;
    static doublereal cr, dr, cw, pl, pr, zk;
    static integer nstack;
    static doublereal wstack[240000], zstack[240000], cl1;
    static integer isteps, jp1;
    static doublereal cw1, cw4, abd, dlc, bod, brd;
    static integer idr, jin;
    static shortint ijstack[240000];

    /* Fortran I/O blocks */
    static cilist io___102 = { 0, 6, 0, 0, 0 };


/*     anisotropic wolff cluster */
/*     yet another time the maximal sizes:(MM) */
/*     stacks are a kind of arrays, seemingly (MM) */
/*     logical prt */
/*     the loop (MM) */
    i__1 = *nsteps;
    for (isteps = 1; isteps <= i__1; ++isteps) {
/*     prt=(isteps.eq.17) */
	ransr_();
/*     choose origin */
	nran_1.nrannr = 1;
	ik = (integer) (nran_1.irn[nran_1.nrannr - 1] * lookup_1.rmul + 
		lookup_1.radd);
	nran_1.nrannr = 2;
	zk = nran_1.irn[nran_1.nrannr - 1] * lookup_1.rmul2 + lookup_1.radd2;
/*     determine sign of cluster (before flip) */
	icsp = mats_1.iz0s[ik - 1];
	jin = 0;
L180:
	if (jin >= mats_1.nin[ik - 1]) {
	    goto L181;
	}
	jp1 = jin + 1;
	if (mats_1.zsbnd[ik + jp1 * 1600 - 1601] > zk) {
	    goto L181;
	}
	jin = jp1;
	icsp = -icsp;
	goto L180;
L181:
	nstack = 0;
	cw1 = 0.;
	nran_1.nrannr = 3;
	brd = lookup_1.rmul3 * nran_1.irn[nran_1.nrannr - 1] + .5;
	brd = -log(brd);
	nran_1.nrannr = 4;
	bod = lookup_1.rmul3 * nran_1.irn[nran_1.nrannr - 1] + .5;
	bod = -log(bod) * lookuc_1.t2p;
	cw4 = cw1 * 4.;
	abd = 0.;
/*     zk,cw4,ik,abd,icsp,brd,bod defined */
L104:
/*     begin determination horiz. bounds of cluster */
	if (mats_1.nin[ik - 1] > 0) {
	    goto L194;
	}
	if (brd >= dats_1.a3) {
	    cl = 0.;
	    cw = dats_1.a3;
	    brd -= dats_1.a3;
	    mats_1.iz0s[ik - 1] = -mats_1.iz0s[ik - 1];
	    goto L193;
	}
	dlc = brd;
	dr = dats_1.a3 - brd;
	++nran_1.nrannr;
	brd = lookup_1.rmul3 * nran_1.irn[nran_1.nrannr - 1] + .5;
	brd = -log(brd);
/*     write(6,889) irn(nrannr),brd */
	if (brd >= dr) {
	    cl = 0.;
	    cw = dats_1.a3;
	    brd -= dr;
	    mats_1.iz0s[ik - 1] = -mats_1.iz0s[ik - 1];
	    goto L193;
	}
	mats_1.nin[ik - 1] = 2;
	cl = zk - dlc;
	if (cl < 0.) {
	    cl += dats_1.a3;
	    mats_1.iz0s[ik - 1] = -mats_1.iz0s[ik - 1];
	    mats_1.zsbnd[ik - 1] = zk + brd;
	    mats_1.zsbnd[ik + 1599] = cl;
	} else {
	    cr = zk + brd;
	    if (cr >= dats_1.a3) {
		mats_1.iz0s[ik - 1] = -mats_1.iz0s[ik - 1];
		mats_1.zsbnd[ik - 1] = cr - dats_1.a3;
		mats_1.zsbnd[ik + 1599] = cl;
	    } else {
		mats_1.zsbnd[ik - 1] = cl;
		mats_1.zsbnd[ik + 1599] = cr;
	    }
	}
	cw = dlc + brd;
/*     calculate new break length */
	++nran_1.nrannr;
	brd = lookup_1.rmul3 * nran_1.irn[nran_1.nrannr - 1] + .5;
	brd = -log(brd);
/*     write(6,889) irn(nrannr),brd */
	goto L193;
L194:
/*     determine size Ising cluster */
	if (jin == 0) {
	    pl = mats_1.zsbnd[ik + mats_1.nin[ik - 1] * 1600 - 1601];
	    dl = zk - pl + dats_1.a3;
	} else {
	    pl = mats_1.zsbnd[ik + jin * 1600 - 1601];
	    dl = zk - pl;
	}
	if (jin == mats_1.nin[ik - 1]) {
	    pr = mats_1.zsbnd[ik - 1];
	    dr = pr + dats_1.a3 - zk;
	} else {
	    pr = mats_1.zsbnd[ik + (jin + 1) * 1600 - 1601];
	    dr = pr - zk;
	}
/*     first find left random cluster reach */
	if (brd >= dl) {
	    cl = pl;
	    brd -= dl;
/*     eliminate left interface */
	    if (jin != 0) {
		i__2 = mats_1.nin[ik - 1];
		for (ii = jin + 1; ii <= i__2; ++ii) {
		    mats_1.zsbnd[ik + (ii - 1) * 1600 - 1601] = mats_1.zsbnd[
			    ik + ii * 1600 - 1601];
/* L183: */
		}
		--jin;
	    } else {
		mats_1.iz0s[ik - 1] = -mats_1.iz0s[ik - 1];
	    }
	    mats_1.nin[ik - 1] = (shortint) (mats_1.nin[ik - 1] - 1);
	} else {
	    dl = brd;
/*     add left interface */
	    cl = zk - brd;
	    mats_1.nin[ik - 1] = (shortint) (mats_1.nin[ik - 1] + 1);
	    if (cl < 0.) {
		cl += dats_1.a3;
		mats_1.zsbnd[ik + mats_1.nin[ik - 1] * 1600 - 1601] = cl;
		mats_1.iz0s[ik - 1] = -mats_1.iz0s[ik - 1];
	    } else {
		++jin;
		i__2 = jin + 1;
		for (ii = mats_1.nin[ik - 1]; ii >= i__2; --ii) {
		    mats_1.zsbnd[ik + ii * 1600 - 1601] = mats_1.zsbnd[ik + (
			    ii - 1) * 1600 - 1601];
/* L184: */
		}
		mats_1.zsbnd[ik + jin * 1600 - 1601] = cl;
	    }
/*     write(6,351) ik,nin(ik),(zsbnd(ik,iii),iii=1,nin(ik)) */
/* 351  format('ad l int: k,nin,zsbnd',2i3,20f6.2) */
	    ++nran_1.nrannr;
	    brd = lookup_1.rmul3 * nran_1.irn[nran_1.nrannr - 1] + .5;
	    brd = -log(brd);
/*     write(6,889) irn(nrannr),brd */
	}
/*  analogous for right hand interval */
	if (brd > dr) {
	    cr = pr;
	    brd -= dr;
	    cw = dl + dr;
/*  eliminate right interface */
	    if (jin != mats_1.nin[ik - 1]) {
		i__2 = mats_1.nin[ik - 1] - 1;
		for (ii = jin + 1; ii <= i__2; ++ii) {
		    mats_1.zsbnd[ik + ii * 1600 - 1601] = mats_1.zsbnd[ik + (
			    ii + 1) * 1600 - 1601];
/* L185: */
		}
	    } else {
		i__2 = mats_1.nin[ik - 1] - 1;
		for (ii = 1; ii <= i__2; ++ii) {
		    mats_1.zsbnd[ik + ii * 1600 - 1601] = mats_1.zsbnd[ik + (
			    ii + 1) * 1600 - 1601];
/* L195: */
		}
		mats_1.iz0s[ik - 1] = -mats_1.iz0s[ik - 1];
	    }
	    mats_1.nin[ik - 1] = (shortint) (mats_1.nin[ik - 1] - 1);
/*     write(6,352) ik,nin(ik),(zsbnd(ik,iii),iii=1,nin(ik)) */
/* 352  format('el r int: k,nin,zsbnd',2i3,20f6.2) */
	} else {
/*  add right interface */
/* so not brd.gt.dr (MM) */
	    dr = brd;
	    cw = dl + dr;
	    cr = zk + brd;
	    mats_1.nin[ik - 1] = (shortint) (mats_1.nin[ik - 1] + 1);
	    if (cr > dats_1.a3) {
		cr -= dats_1.a3;
		for (ii = mats_1.nin[ik - 1]; ii >= 2; --ii) {
		    mats_1.zsbnd[ik + ii * 1600 - 1601] = mats_1.zsbnd[ik + (
			    ii - 1) * 1600 - 1601];
/* L187: */
		}
		mats_1.zsbnd[ik - 1] = cr;
		mats_1.iz0s[ik - 1] = -mats_1.iz0s[ik - 1];
	    } else {
		i__2 = jin + 2;
		for (ii = mats_1.nin[ik - 1]; ii >= i__2; --ii) {
		    mats_1.zsbnd[ik + ii * 1600 - 1601] = mats_1.zsbnd[ik + (
			    ii - 1) * 1600 - 1601];
/* L186: */
		}
		mats_1.zsbnd[ik + (jin + 1) * 1600 - 1601] = cr;
	    }
	    ++nran_1.nrannr;
	    brd = lookup_1.rmul3 * nran_1.irn[nran_1.nrannr - 1] + .5;
	    brd = -log(brd);
/*     write(6,889) irn(nrannr),brd */
/* 889  format('breakd.:'i12,2f20.14) */
	}
/*  end determination horiz. bounds of flipped part */
/*  write in stack */
L193:
	++nstack;
	ijstack[nstack - 1] = (shortint) ik;
	zstack[nstack - 1] = cl;
	wstack[nstack - 1] = cw;
/*     write(6,888) nstack,ik,cl,cw */
/* 888  format('stackwr:'2i6,2f20.14) */
L192:
	abd += bod;
	if (abd >= cw4) {
	    bod = abd - cw4;
	    goto L161;
	}
	idr = (integer) (abd / cw1);
	zk = cl1 + abd - idr * cw1;
	if (zk > dats_1.a3) {
	    zk -= dats_1.a3;
	}
	ik = ngcb_1.ngc[idr + 1 + (ij << 2) - 5];
/*     find sign of nb */
	isnb = mats_1.iz0s[ik - 1];
	jin = 0;
/*     while-like construction(MM) */
L188:
	if (jin >= mats_1.nin[ik - 1]) {
	    goto L189;
	}
	jp1 = jin + 1;
	if (mats_1.zsbnd[ik + jp1 * 1600 - 1601] > zk) {
	    goto L189;
	}
	jin = jp1;
	isnb = -isnb;
	goto L188;
L189:
/* end find sign */
/* first determine interval to next perp bond */
	++nran_1.nrannr;
	bod = lookup_1.rmul3 * nran_1.irn[nran_1.nrannr - 1] + .5;
	bod = -log(bod) * lookuc_1.t2p;
	if (isnb == icsp) {
	    goto L104;
	}
/*     104: just before determining the horizontal bound (MM) */
	goto L192;
/*     192: just after write in stack (MM) */
L161:
	if (nstack == 0) {
	    goto L105;
	}
/* read from stack */
	ij = ijstack[nstack - 1];
	cl1 = zstack[nstack - 1];
	cw1 = wstack[nstack - 1];
	cw4 = cw1 * 4.;
	--nstack;
/*     write(6,981) isteps,nstack,ij,cl1,cw1 */
/*     981  format ('str isteps,nstack ij cl cw',3i4,2f8.3) */
	abd = 0.;
	goto L192;
L105:
/* stack was empty (MM) */
	if (nran_1.nrannr > lookuc_1.maxr) {
	    lookuc_1.maxr = nran_1.nrannr;
	    if (nran_1.nrannr > 480000) {
		s_wsle(&io___102);
		do_lio(&c__9, &c__1, " **** nrannr>mxrn:", 18L);
		do_lio(&c__3, &c__1, (char *)&nran_1.nrannr, (ftnlen)sizeof(
			integer));
		do_lio(&c__3, &c__1, (char *)&c_b98, (ftnlen)sizeof(integer));
		e_wsle();
/* L2: */
		s_stop("", 0L);
	    }
	}
/* L150: */
    }
/*     150 end of main loop of this subroutine (MM) */
    return 0;
} /* mcw_ */

/*     end of subroutine mcw */
/*  zk intoduced everywhere; OK for both par and perp bonds? */
/* parallel parms: append 1 */
/* vermijd interf parms in par and perp sectie */
/* interfaces (grain boundaries) everywhere updated? */
/* check precise bounds of intervals generated by bod and brd */
/* --------------------------------------------------------------- */
/* core */


/* Subroutine */ int core_(sn, ncor)
doublereal *sn;
integer *ncor;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static doublereal z__, s2;
    static integer ij;
    static doublereal st, vz;
    static integer ij1, ij2, ni1, in1, in2, ni2, jn2, jn1, inb, nij;
    static doublereal amt, bsn;
    static integer nit, int__, izs;
    static doublereal sps;
    static integer isp2;

/*     version for cylinder ising model */
/*     The maximal sizes...(MM) */
/*     the neighbors: (MM) */
    /* Parameter adjustments */
    --sn;

    /* Function Body */
    *ncor = 9;
/*     if (ncor.eq.9) goto 10 */
/*     sample correlation function over half the system size */
/* L10: */
/*     sample magnetization and nr of interfaces */
    nit = 0;
    amt = 0.;
/*     n12 is size of system (MM) */
    i__1 = dats_1.n12;
    for (ij = 1; ij <= i__1; ++ij) {
	nij = mats_1.nin[ij - 1];
	if (nij > lookuc_1.maxi) {
	    lookuc_1.maxi = nij;
	}
	nit += nij;
	izs = mats_1.iz0s[ij - 1];
	amt += dats_1.a3 * izs;
	if (nij > 0) {
	    vz = mats_1.zsbnd[ij + nij * 1600 - 1601];
	}
	i__2 = nij;
	for (int__ = 1; int__ <= i__2; ++int__) {
	    z__ = mats_1.zsbnd[ij + int__ * 1600 - 1601];
	    amt += (z__ - vz) * izs;
	    izs = -izs;
	    vz = z__;
/* L19: */
	}
/* L18: */
    }
/*     another loop of the whole system (MM) */
    sps = 0.;
    i__1 = dats_1.n12;
    for (ij1 = 1; ij1 <= i__1; ++ij1) {
	ni1 = mats_1.nin[ij1 - 1];
	for (inb = 1; inb <= 2; ++inb) {
	    in1 = 0;
	    in2 = 0;
	    ij2 = ngcb_1.ngc[inb + (ij1 << 2) - 5];
	    ni2 = mats_1.nin[ij2 - 1];
	    isp2 = mats_1.iz0s[ij1 - 1] * mats_1.iz0s[ij2 - 1];
	    sps += dats_1.a3 * isp2;
	    isp2 += isp2;
	    if (0 == ni1) {
		goto L21;
	    }
	    if (0 == ni2) {
		goto L22;
	    }
L20:
	    if (mats_1.zsbnd[ij1 + (in1 + 1) * 1600 - 1601] > mats_1.zsbnd[
		    ij2 + (in2 + 1) * 1600 - 1601]) {
		++in2;
		sps += isp2 * mats_1.zsbnd[ij2 + in2 * 1600 - 1601];
		isp2 = -isp2;
		if (in2 == ni2) {
		    goto L22;
		}
	    } else {
		++in1;
		sps += isp2 * mats_1.zsbnd[ij1 + in1 * 1600 - 1601];
		isp2 = -isp2;
		if (in1 == ni1) {
		    goto L21;
		}
	    }
	    goto L20;
L21:
	    i__2 = ni2;
	    for (jn2 = in2 + 1; jn2 <= i__2; ++jn2) {
		sps += isp2 * mats_1.zsbnd[ij2 + jn2 * 1600 - 1601];
		isp2 = -isp2;
/* L23: */
	    }
	    goto L25;
L22:
	    i__2 = ni1;
	    for (jn1 = in1 + 1; jn1 <= i__2; ++jn1) {
		sps += isp2 * mats_1.zsbnd[ij1 + jn1 * 1600 - 1601];
		isp2 = -isp2;
/* L24: */
	    }
L25:
/* L26: */
	    ;
	}
/* L27: */
    }
/*     I've added this 27, both loops pointed to 26, */
/*     which of course is possible but emacs doesn't like it (MM) */
/*     the things this program calculates: (MM) */
    bsn = nit * dats_1.wn;
    sps *= dats_1.wn;
    st = amt * dats_1.wn;
    s2 = st * st;
    sn[1] = s2;
    sn[2] = s2 * s2;
    sn[3] = bsn;
    sn[4] = bsn * bsn;
    sn[5] = sn[3] * sn[4];
    sn[6] = sn[4] * sn[4];
    sn[7] = sn[1] * sn[3];
    sn[8] = sn[2] * sn[3];
    sn[9] = sps;
    return 0;
} /* core_ */

/*     end of core */
/* -------------------------------------------------------------- */
/* setpqr */
/* this subrouting probably sets p q (to zero?) (MM) */
/* Subroutine */ int setpqr_()
{
    static integer k, l;

/*     npmax?(MM) */
/*     the maximal sizes (MM) */
    for (k = 1; k <= 9; ++k) {
	pqrh_1.p[k - 1] = (float)0.;
	for (l = 1; l <= 9; ++l) {
	    pqrh_1.q[k + l * 9 - 10] = (float)0.;
/* L1: */
	}
/* L2: */
    }
/*     I've added '2' (MM) */
    for (k = 1; k <= 2000; ++k) {
	pqrh_1.ipr[k - 1] = 0;
/* L4: */
    }
/*     I've added this continue (MM) */
    return 0;
} /* setpqr_ */

/*     end of setpqr */
/* ---------------------------------------------------------- */
/* addpqr */
/* p,q are things which are greatened every iterations, by means */
/* of this funcion. (MM) */
/* Subroutine */ int addpqr_(n1, n2, ntop)
integer *n1, *n2, *ntop;
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int core_();
    static integer j, k;
    static doublereal sn[9];

/*     again.. (MM) */
    core_(sn, ntop);
    i__1 = *ntop;
    for (j = 1; j <= i__1; ++j) {
	pqrh_1.p[j - 1] += sn[j - 1];
	i__2 = *ntop;
	for (k = 1; k <= i__2; ++k) {
	    pqrh_1.q[j + k * 9 - 10] += sn[j - 1] * sn[k - 1];
/* L1: */
	}
/* L2: */
    }
/*     i've added 2 (MM) */
    return 0;
} /* addpqr_ */

/*     end of addpqr */
/* -------------------------------------------------------------- */
/* part */

/* Subroutine */ int part_(ip, nc, ntop)
integer *ip, *nc, *ntop;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, nv;

/* in part: */
    i__1 = *ntop;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pqrh_1.pp[i__ + *ip * 9 - 10] = pqrh_1.p[i__ - 1];
	pqrh_1.p[i__ - 1] = 0.;
/* L100: */
    }
    if (*ip == 1) {
	nv = 0;
    }
    pqrh_1.ipr[*ip - 1] = *nc - nv;
    nv = *nc;
/* end */
    return 0;
} /* part_ */

/*     end of part */
/* ----------------------------------------------------- */
/* pqr */


/* Subroutine */ int pqr_(ncup, ncount, npart, n1, n2, a3)
integer *ncup, *ncount, *npart, *n1, *n2;
doublereal *a3;
{
    /* Initialized data */

    static integer ndiv[12] = { 2,5,10,20,50,100,200,500,1000,2000,5000,10000 
	    };

    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Local variables */
    extern /* Subroutine */ int avdv_();
    static integer ntog, i__, j, k, npdiv, nrdiv;
    static doublereal rp[9000]	/* was [9][1000] */;
    static integer ind, npd, jpr[1000];
    static logical prt;

    prt = FALSE_;
    i__1 = *ncup;
    for (j = 1; j <= i__1; ++j) {
	pqrh_1.p[j - 1] = 0.;
	i__2 = *npart;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    pqrh_1.p[j - 1] += pqrh_1.pp[j + i__ * 9 - 10];
/* L11: */
	}
	pqrh_1.p[j - 1] /= *ncount;
/* L5: */
    }
    i__1 = *ncup;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *ncup;
	for (k = 1; k <= i__2; ++k) {
	    pqrh_1.q[j + k * 9 - 10] /= *ncount;
	    cbh_1.b[j + k * 9 - 10] = pqrh_1.q[j + k * 9 - 10] - pqrh_1.p[j - 
		    1] * pqrh_1.p[k - 1];
/* L1: */
	}
/* L2: */
    }
/* (MM) */
    i__1 = *npart;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (pqrh_1.ipr[i__ - 1] == 0) {
	    goto L20;
	}
	i__2 = *ncup;
	for (j = 1; j <= i__2; ++j) {
	    pqrh_1.pp[j + i__ * 9 - 10] /= pqrh_1.ipr[i__ - 1];
/* L21: */
	}
/* (MM) */
L20:
	;
    }
    avdv_(npart, ncup, pqrh_1.ipr, pqrh_1.pp, pqrh_1.p, &prt, n1, n2, a3);
    nrdiv = 1;
L30:
    npdiv = *npart / ndiv[nrdiv - 1];
    prt = nrdiv == 1;
    if (npdiv <= 3) {
	return 0;
    }
/*     this is the only way to exit this subroutine */
/*     the following is done as long as necessary (MM) */
    i__1 = npdiv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	jpr[i__ - 1] = 0;
	i__2 = *ncup;
	for (j = 1; j <= i__2; ++j) {
	    rp[j + i__ * 9 - 10] = 0.;
/* L31: */
	}
/* L29: */
    }
/*  (MM) */
    npd = 0;
    ind = 0;
    i__1 = npdiv;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ntog = ndiv[nrdiv - 1];
	i__2 = ntog;
	for (k = 1; k <= i__2; ++k) {
	    ++ind;
	    jpr[i__ - 1] += pqrh_1.ipr[ind - 1];
	    i__3 = *ncup;
	    for (j = 1; j <= i__3; ++j) {
		rp[j + i__ * 9 - 10] += pqrh_1.pp[j + ind * 9 - 10] * 
			pqrh_1.ipr[ind - 1];
/* L34: */
	    }
/* L33: */
	}
	if (jpr[i__ - 1] == 0) {
	    goto L32;
	}
	++npd;
	i__2 = *ncup;
	for (j = 1; j <= i__2; ++j) {
	    rp[j + i__ * 9 - 10] /= jpr[i__ - 1];
/* L35: */
	}
L32:
	;
    }
    avdv_(&npdiv, ncup, jpr, rp, pqrh_1.p, &prt, n1, n2, a3);
    ++nrdiv;
    goto L30;
} /* pqr_ */

/*     end pqr */
/* -------------------------------------------------------- */
/* avdv */

/* Subroutine */ int avdv_(npart, ncup, ipr, pp, p, prt, n1, n2, a3)
integer *npart, *ncup, *ipr;
doublereal *pp, *p;
logical *prt;
integer *n1, *n2;
doublereal *a3;
{
    /* Format strings */
    static char fmt_70[] = "(\0020npart=\002,i8,\002 : partial results suppr\
essed\002//)";
    static char fmt_60[] = "(\0020npart=\002,i8,\002 partial results for ave\
rages\002//)";
    static char fmt_61[] = "(\002 sum terms\002,9(i8,4x))";
    static char fmt_62[] = "(//\002 \002)";
    static char fmt_64[] = "(\002 \002,i3,i8,10f12.6)";
    static char fmt_82[] = "(/\002  nc     average    st. dev.    corr. c\
.\002)";
    static char fmt_80[] = "(\002 \002,i3,3f14.8,i3,3f12.6,3x,3f12.6)";

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();
    double sqrt();

    /* Local variables */
    static doublereal aver[20], sdev[9], vdev, corr[9];
    static integer nctr, i__, j, k;
    static doublereal b3, b4, sdevj1, corrj1, sc;
    static integer np;
    static doublereal sk;
    static integer n1h, n2h;
    static doublereal dev, vcu;

    /* Fortran I/O blocks */
    static cilist io___148 = { 0, 6, 0, fmt_70, 0 };
    static cilist io___149 = { 0, 6, 0, fmt_60, 0 };
    static cilist io___150 = { 0, 6, 0, fmt_61, 0 };
    static cilist io___152 = { 0, 6, 0, fmt_62, 0 };
    static cilist io___154 = { 0, 6, 0, fmt_64, 0 };
    static cilist io___163 = { 0, 6, 0, fmt_82, 0 };
    static cilist io___164 = { 0, 6, 0, fmt_80, 0 };
    static cilist io___165 = { 0, 8, 0, fmt_80, 0 };
    static cilist io___171 = { 0, 6, 0, fmt_80, 0 };
    static cilist io___172 = { 0, 8, 0, fmt_80, 0 };


/*     real*8 pp(mxcr,npmax),p(mxcr),sdev(mxcr) */
/*     double precision wn,dev,vdev,sk,sc */
    /* Parameter adjustments */
    --p;
    pp -= 10;
    --ipr;

    /* Function Body */
    n1h = *n1 / 2;
    n2h = *n2 / 2;
    nctr = *ncup;
    vcu = *n1 * *n2 * *a3;
    if (nctr > 6) {
	nctr = 6;
    }
    if (*npart >= 100) {
	s_wsfe(&io___148);
	do_fio(&c__1, (char *)&(*npart), (ftnlen)sizeof(integer));
	e_wsfe();
	goto L71;
    }
    s_wsfe(&io___149);
    do_fio(&c__1, (char *)&(*npart), (ftnlen)sizeof(integer));
    e_wsfe();
    s_wsfe(&io___150);
    i__1 = nctr;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_fio(&c__1, (char *)&i__, (ftnlen)sizeof(integer));
    }
    e_wsfe();
    s_wsfe(&io___152);
    e_wsfe();
    i__1 = *npart;
    for (k = 1; k <= i__1; ++k) {
	s_wsfe(&io___154);
	do_fio(&c__1, (char *)&k, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&ipr[k], (ftnlen)sizeof(integer));
	i__2 = nctr;
	for (j = 1; j <= i__2; ++j) {
	    do_fio(&c__1, (char *)&pp[j + k * 9], (ftnlen)sizeof(doublereal));
	}
	e_wsfe();
/* L63: */
    }
    if (*npart <= 3) {
	return 0;
    }
L71:
    i__1 = *ncup;
    for (j = 1; j <= i__1; ++j) {
	np = 0;
	sk = 0.;
	sc = 0.;
	dev = 0.;
	i__2 = *npart;
	for (k = 1; k <= i__2; ++k) {
	    if (ipr[k] == 0) {
		goto L12;
	    }
/* break out of this do-loop (MM) */
	    vdev = dev;
	    dev = pp[j + k * 9] - p[j];
	    sk += dev * dev;
	    sc += dev * vdev;
	    ++np;
/* L11: */
	}
L12:
	if (np < 3) {
	    return 0;
	}
	sdev[j - 1] = sqrt(sk / (np * (np - 1)));
	if (sk != 0.) {
	    corr[j - 1] = np * sc / (sk * (np - 1));
	} else {
	    corr[j - 1] = 0.;
	}
/* L10: */
    }
    s_wsfe(&io___163);
    e_wsfe();
    i__1 = *ncup;
    for (j = 1; j <= i__1; ++j) {
	s_wsfe(&io___164);
	do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&p[j], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&sdev[j - 1], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&corr[j - 1], (ftnlen)sizeof(doublereal));
	e_wsfe();
	if (*prt) {
	    s_wsfe(&io___165);
	    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&p[j], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&sdev[j - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&corr[j - 1], (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
/* L15: */
    }
    b3 = vcu;
    b4 = b3 * vcu;
/* Computing 2nd power */
    d__1 = p[1];
    aver[0] = d__1 * d__1 / p[2];
/* q */
/* Computing 2nd power */
    d__1 = p[3];
    aver[1] = vcu * (p[4] - d__1 * d__1);
/* c */
/* Computing 3rd power */
    d__1 = p[3], d__2 = d__1;
    aver[2] = b3 * (p[5] - p[4] * 3. * p[3] + d__2 * (d__1 * d__1) * 2.);
/* c' */
/* Computing 2nd power */
    d__1 = p[3];
/* Computing 2nd power */
    d__2 = p[4];
/* Computing 4th power */
    d__3 = p[3], d__3 *= d__3;
    aver[3] = b4 * (p[6] - p[5] * 4. * p[3] + p[4] * 12. * (d__1 * d__1) - 
	    d__2 * d__2 * 3. - d__3 * d__3 * 6.);
/* c'' */
    aver[4] = p[7] - p[3] * p[1];
/* em2 */
    aver[5] = p[8] - p[3] * p[2];
/* em4 */
    aver[6] = p[7] * 2. / p[1] - p[8] / p[2] - p[3];
/* qp */
    for (i__ = 1; i__ <= 7; ++i__) {
	j = *ncup + i__;
	np = 0;
	sk = 0.;
	sc = 0.;
	dev = 0.;
	i__1 = *npart;
	for (k = 1; k <= i__1; ++k) {
	    vdev = dev;
	    switch ((int)i__) {
		case 1:  goto L30;
		case 2:  goto L31;
		case 3:  goto L32;
		case 4:  goto L33;
		case 5:  goto L34;
		case 6:  goto L35;
		case 7:  goto L36;
	    }
L30:
	    if (pp[k * 9 + 1] == 0.) {
		dev = .33333333333333331 - aver[0];
	    } else {
/* Computing 2nd power */
		d__1 = pp[k * 9 + 1];
		dev = d__1 * d__1 / pp[k * 9 + 2] - aver[0];
	    }
	    goto L50;
L31:
/* Computing 2nd power */
	    d__1 = pp[k * 9 + 3];
	    dev = vcu * (pp[k * 9 + 4] - d__1 * d__1) - aver[1];
	    goto L50;
L32:
/* Computing 3rd power */
	    d__1 = pp[k * 9 + 3], d__2 = d__1;
	    dev = b3 * (pp[k * 9 + 5] - pp[k * 9 + 4] * 3. * pp[k * 9 + 3] + 
		    d__2 * (d__1 * d__1) * 2.) - aver[2];
	    goto L50;
L33:
/* Computing 2nd power */
	    d__1 = pp[k * 9 + 3];
/* Computing 2nd power */
	    d__2 = pp[k * 9 + 4];
/* Computing 4th power */
	    d__3 = pp[k * 9 + 3], d__3 *= d__3;
	    dev = b4 * (pp[k * 9 + 6] - pp[k * 9 + 5] * 4. * pp[k * 9 + 3] + 
		    pp[k * 9 + 4] * 12. * (d__1 * d__1) - d__2 * d__2 * 3. - 
		    d__3 * d__3 * 6.) - aver[3];
	    goto L50;
L34:
	    dev = pp[k * 9 + 7] - pp[k * 9 + 3] * pp[k * 9 + 1] - aver[4];
	    goto L50;
L35:
	    dev = pp[k * 9 + 8] - pp[k * 9 + 3] * pp[k * 9 + 2] - aver[5];
	    goto L50;
L36:
	    if (pp[k * 9 + 1] == 0.) {
		dev = -aver[6];
	    } else {
		dev = pp[k * 9 + 7] * 2 / pp[k * 9 + 1] - pp[k * 9 + 8] / pp[
			k * 9 + 2] - pp[k * 9 + 3] - aver[6];
	    }
L50:
	    sk += dev * dev;
	    sc += dev * vdev;
	    ++np;
/* L21: */
	}
	sdevj1 = sqrt(sk / (np * (np - 1)));
	if (sk != 0.) {
	    corrj1 = np * sc / (sk * (np - 1));
	} else {
	    corrj1 = 0.;
	}
	s_wsfe(&io___171);
	do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&aver[i__ - 1], (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&sdevj1, (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&corrj1, (ftnlen)sizeof(doublereal));
	e_wsfe();
	if (*prt) {
	    s_wsfe(&io___172);
	    do_fio(&c__1, (char *)&j, (ftnlen)sizeof(integer));
	    do_fio(&c__1, (char *)&aver[i__ - 1], (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&sdevj1, (ftnlen)sizeof(doublereal));
	    do_fio(&c__1, (char *)&corrj1, (ftnlen)sizeof(doublereal));
	    e_wsfe();
	}
/* L22: */
    }
/* L20: */
/*     write(6,81) aver(1),aver(2),aver(3),aver(4) */
/* 81   format(' q=',f10.6,' c=',f10.5,' em2=',f12.7,' em4=',f12.7) */
/* where is this one for ? (MM)? */
    return 0;
} /* avdv_ */

/*     end avdv (MM) */
/* --------------------------------------------------------------- */
/* ransi */
/* random numbers (MM) */
/* Subroutine */ int ransi_0_(n__, iseed)
int n__;
integer *iseed;
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer ipnf1, ipnf2, ipnt1, ipnt2, i__, k, l, k1;

/* sequential version */
/*     shift register random generator with very long period */
    switch(n__) {
	case 1: goto L_ransr;
	}

    k = (*iseed << 1) + 387420489;
    k1 = *iseed * 1313131;
    k1 -= k1 / 2796203 * 2796203;
    for (i__ = 1; i__ <= 9689; ++i__) {
	k *= 32781;
	k1 *= 125;
	k1 -= k1 / 2796203 * 2796203;
	ransrb_1.ir1[i__ - 1] = k + k1 * 8193;
/* L100: */
    }
    for (i__ = 1; i__ <= 127; ++i__) {
	k *= 32781;
	k1 *= 125;
	k1 -= k1 / 2796203 * 2796203;
	ransrb_1.ir2[i__ - 1] = k + k1 * 4099;
/* L102: */
    }
    for (i__ = 1; i__ <= 9689; ++i__) {
	ransrb_1.inxt1[i__ - 1] = (shortint) (i__ + 1);
/* L101: */
    }
/* MM */
    ransrb_1.inxt1[9688] = 1;
    ipnt1 = 1;
    ipnf1 = 472;
    for (i__ = 1; i__ <= 127; ++i__) {
	ransrb_1.inxt2[i__ - 1] = (shortint) (i__ + 1);
/* L103: */
    }
/* MM */
    ransrb_1.inxt2[126] = 1;
    ipnt2 = 1;
    ipnf2 = 31;
    return 0;
/* -----------------entry ransr-------------------------------- */

L_ransr:
/*     calculate n random numbers */
    i__1 = nran_2.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	l = ransrb_1.ir1[ipnt1 - 1] ^ ransrb_1.ir1[ipnf1 - 1];
	k = ransrb_1.ir2[ipnt2 - 1] ^ ransrb_1.ir2[ipnf2 - 1];
	nran_2.irn[i__ - 1] = k ^ l;
	ransrb_1.ir1[ipnt1 - 1] = l;
	ipnt1 = ransrb_1.inxt1[ipnt1 - 1];
	ipnf1 = ransrb_1.inxt1[ipnf1 - 1];
	ransrb_1.ir2[ipnt2 - 1] = k;
	ipnt2 = ransrb_1.inxt2[ipnt2 - 1];
	ipnf2 = ransrb_1.inxt2[ipnf2 - 1];
/* L200: */
    }
    return 0;
} /* ransi_ */

/* Subroutine */ int ransi_(iseed)
integer *iseed;
{
    return ransi_0_(0, iseed);
    }

/* Subroutine */ int ransr_()
{
    return ransi_0_(1, (integer *)0);
    }

/*     end ransi/ransr (MM) */
/* ----------------------------------------------------------------- */
/* ctijd */

doublereal ctijd_(a)
real *a;
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    static doublereal t;
    extern doublereal etime_();
    static real ta[2];

/*     function ctijd returns time in seconds; a is a dummy argument */
    t = etime_(ta);
/*     ctijd = ta(1) */
    ret_val = t;
    return ret_val;
} /* ctijd_ */

/* Main program alias */ int anlim_ () { MAIN__ (); }
