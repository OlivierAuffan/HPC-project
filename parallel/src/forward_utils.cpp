#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "forward.h"
#include "shalw.h"

double hFil_forward(int t, int i, int j)
{
    if (t <= 2)
	return HPHY(t, i, j);
    return HPHY(t - 1, i, j) +
	alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
}

double hFil_forward_simd(int t, int i, int j)
{
    return HPHY(t - 1, i, j) +
	alpha * (HFIL(t - 1, i, j) - 2 * HPHY(t - 1, i, j) + HPHY(t, i, j));
}

/////////////////////////////////////////////////////////////////////////////////

double uFil_forward(int t, int i, int j)
{
    if (t <= 2)
	return UPHY(t, i, j);
    return UPHY(t - 1, i, j) +
	alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j));
}

double uFil_forward_simd(int t, int i, int j)
{
    return UPHY(t - 1, i, j) +
	alpha * (UFIL(t - 1, i, j) - 2 * UPHY(t - 1, i, j) + UPHY(t, i, j));
}

/////////////////////////////////////////////////////////////////////////////////

double vFil_forward(int t, int i, int j)
{
    if (t <= 2)
	return VPHY(t, i, j);
    return VPHY(t - 1, i, j) +
	alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j));
}

double vFil_forward_simd(int t, int i, int j)
{
    return VPHY(t - 1, i, j) +
	alpha * (VFIL(t - 1, i, j) - 2 * VPHY(t - 1, i, j) + VPHY(t, i, j));
}

/////////////////////////////////////////////////////////////////////////////////

double hPhy_forward(int t, int i, int j)
{
    double c, d;

    c = 0.;
    if (i > 0)
	c = UPHY(t - 1, i - 1, j);

    d = 0.;
    if (j < size_y - 1)
	d = VPHY(t - 1, i, j + 1);

    return HFIL(t - 1, i, j) -
	dt * hmoy *
	((UPHY(t - 1, i, j) - c) / dx + (d - VPHY(t - 1, i, j)) / dy);
}

double hPhy_forward_simd(int t, int i, int j)
{
    double c, d;

    c = UPHY(t - 1, i - 1, j);
    d = VPHY(t - 1, i, j + 1);

    return HFIL(t - 1, i, j) -
	dt * hmoy *
	((UPHY(t - 1, i, j) - c) / dx + (d - VPHY(t - 1, i, j)) / dy);
}

/////////////////////////////////////////////////////////////////////////////////

double uPhy_forward(int t, int i, int j)
{
    double b, e, f, g;

    if (i == size_x - 1)
	return 0.;

    b = 0.;
    if (i < size_x - 1)
	b = HPHY(t - 1, i + 1, j);

    e = 0.;
    if (j < size_y - 1)
	e = VPHY(t - 1, i, j + 1);

    f = 0.;
    if (i < size_x - 1)
	f = VPHY(t - 1, i + 1, j);

    g = 0.;
    if (i < size_x - 1 && j < size_y - 1)
	g = VPHY(t - 1, i + 1, j + 1);

    return UFIL(t - 1, i, j) +
	    dt * ((-grav / dx) * (b - HPHY(t - 1, i, j)) +
		  (pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g) -
		  (dissip * UFIL(t - 1, i, j)));
}

double uPhy_forward_simd(int t, int i, int j)
{
    double b, e, f, g;

    b = HPHY(t - 1, i + 1, j);
    e = VPHY(t - 1, i, j + 1);
    f = VPHY(t - 1, i + 1, j);
    g = VPHY(t - 1, i + 1, j + 1);

    return UFIL(t - 1, i, j) +
	dt * ((-grav / dx) * (b - HPHY(t - 1, i, j)) +
	      (pcor / 4.) * (VPHY(t - 1, i, j) + e + f + g) -
	      (dissip * UFIL(t - 1, i, j)));
}

/////////////////////////////////////////////////////////////////////////////////

double vPhy_forward(int t, int i, int j)
{
    double c, d, e, f;

    if (j == 0)
	return 0.;

    c = 0.;
    if (j > 0)
	c = HPHY(t - 1, i, j - 1);

    d = 0.;
    if (i > 0 && j > 0)
	d = UPHY(t - 1, i - 1, j - 1);

    e = 0.;
    if (i > 0)
	e = UPHY(t - 1, i - 1, j);

    f = 0.;
    if (j > 0)
	f = UPHY(t - 1, i, j - 1);

    return VFIL(t - 1, i, j) +
	dt * ((-grav / dy) * (HPHY(t - 1, i, j) - c) -
	      (pcor / 4.) * (d + e + f + UPHY(t - 1, i, j)) -
	      (dissip * VFIL(t - 1, i, j)));
}

double vPhy_forward_simd(int t, int i, int j)
{
    double c, d, e, f;

    c = HPHY(t - 1, i, j - 1);
    d = UPHY(t - 1, i - 1, j - 1);
    e = UPHY(t - 1, i - 1, j);
    f = UPHY(t - 1, i, j - 1);

    return VFIL(t - 1, i, j) +
	dt * ((-grav / dy) * (HPHY(t - 1, i, j) - c) -
	      (pcor / 4.) * (d + e + f + UPHY(t - 1, i, j)) -
	      (dissip * VFIL(t - 1, i, j)));
}

