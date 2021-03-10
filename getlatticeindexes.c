/*****************************************************************************/
/*                                                                           */
/* serpent 2 (beta-version) : getlatticeindexes.c                            */
/*                                                                           */
/* Created:       2010/10/08 (JLe)                                           */
/* Last modified: 2016/10/03 (JLe)                                           */
/* Version:       2.1.28                                                     */
/*                                                                           */
/* Description: Returns indexes in squre and hexagonal lattices              */
/*                                                                           */
/* Comments: - Taken from Serpent 1.1.0                                      */
/*                                                                           */
/*****************************************************************************/

#include "header.h"
#include "locations.h"

#define FUNCTION_NAME "GetLatticeIndexes:"

/*****************************************************************************/

long GetLatticeIndexes(double px, double py, double pz, double x0, double y0,
                       double z0, long *tri, long *i, long *j, long *k, long type)
{
  double x, y, n1, n2, n3;
  double mid1, mid2, mid3;

  if ((type == LAT_TYPE_S) || (type == LAT_TYPE_INFS))
    {
      /* Neliöhila, triviaalitapaus */

      *i = (long)rint(x0/px);
      *j = (long)rint(y0/py);
      *k = (long)rint(z0/pz);

      return 0;
    }
  else if ((type == LAT_TYPE_HXT) || (type == LAT_TYPE_HYT))
    {
      if (type == LAT_TYPE_HYT)
        {
          y = x0/px;
          x = y0/py;
        }
      else
        {
          x = x0/px;
          y = y0/py;
        }

      n1 =        2.0*x - 0.5;
      n2 = -x + SQRT3*y - 0.5;
      n3 = -x - SQRT3*y - 0.5;

      mid1 = rint(n1);
      mid2 = rint(n2);
      mid3 = rint(n3);

      *i = (long)floor(0.5 + (mid1 - mid2)/3.0);
      *j = (long)floor(0.5 + (mid2 - mid3)/3.0);
      *k = (long)rint(z0/pz);

	  n1 = floor(0.5 + (mid1 - mid2)/3.0);
	  n2 = floor(0.5 + (mid2 - mid3)/3.0);

	  if (type == LAT_TYPE_HXT)
      {
        x = x0 - (n1 + 0.5 * n2)*px;
        y = y0 - (0.5*SQRT3)*py*n2;

        if ( (y/x) > (1/SQRT3) )
		{
		if ( x > 0) {*tri = 0;}
		else {*tri = 3;}
		}
	  else if ( (y/x) < (-1/SQRT3) )
		{
		if (x > 0) {*tri = 2;}
        else {*tri = 5;}
		}
	  else
		{
		if ( x > 0) {*tri=1;}
		else {*tri=4;}
		}

      }
      else
      {
        x = x0 - (0.5*SQRT3)*px*n2;
        y = y0 - (n1 + 0.5 * n2)*px;

      if (((y/x) > SQRT3) || ((y/x) < -1*SQRT3 ))
		{
		if ( y > 0) {*tri = 1;}
		else {*tri = 4;}
		}
	  else
		{
		if ( (y/x) > 0)
          {
          if (y > 0) {*tri = 0;}
          else {*tri = 3;}
		  }
		else
		  {
		  if (y > 0) {*tri = 2;}
		  else {*tri = 5;}
		  }
		}
      }
      return 0;
    }
  else
    {
      /* Kolmiohila, käytetään R. Mattilan johtamaa kaavaa */

      if ((type == LAT_TYPE_HY) || (type == LAT_TYPE_INFHY))
        {
          /* Y-tyyppi -> vaihdetaan koordinaatit */

          y = x0/px;
          x = y0/py;
        }
      else
        {
          /* X-tyyppi -> käytetään suoraan */

          x = x0/px;
          y = y0/py;
        }

      n1 =        2.0*x - 0.5;
      n2 = -x + SQRT3*y - 0.5;
      n3 = -x - SQRT3*y - 0.5;

      mid1 = rint(n1);
      mid2 = rint(n2);
      mid3 = rint(n3);

      *i = (long)floor(0.5 + (mid1 - mid2)/3.0);
      *j = (long)floor(0.5 + (mid2 - mid3)/3.0);

      /* Z-indeksi */

      *k = (long)rint(z0/pz);

      return 0;
    }

  fprintf(errp, "%s Invalid lattice type %ld.\n", FUNCTION_NAME, type);

  exit(-1);
}

/*****************************************************************************/
