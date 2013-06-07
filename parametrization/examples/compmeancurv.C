/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of GoTools.
 *
 * GoTools is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * GoTools is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with GoTools. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using GoTools.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the GoTools library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#include <iostream>
#include <vector>

using namespace std;


// Taking second order differences
void comp_finite_diff(int M, int N, int dir,
		      const vector<double>& pts,
		      vector<double>& diffs)
{
    double h;
    int off;
    if (dir == 0) {
	// U-direction
	h = 1.0/(M-1.0);
	off = 1;
    } else {
	// V-direction
	h = 1.0/(N-1.0);
	off = M;
    }

    diffs.clear();
    diffs.resize(M*N, 0.0);
    for (int i = 1; i < M-1; ++i) {
	for (int j = 1; j < N-1; ++j) {
	    diffs[i*M + j] = pts[i*M + j + off] - pts[i*M + j - off];
	    diffs[i*M + j] /= 2*h;
	}
    }
    
}


int main()
{
    // Read the grid
    int M, N;
    cin >> M >> N;
    vector<double> f(M*N);
    for (int i = 0; i < M*N; ++i) {
	cin >> f[i];
    }

    vector<double> fu, fv, fuu, fuv, fvv, sum;

    comp_finite_diff(M, N, 0, f, fu);
    comp_finite_diff(M, N, 1, f, fv);
    comp_finite_diff(M, N, 0, fu, fuu);
    comp_finite_diff(M, N, 1, fu, fuv);
    comp_finite_diff(M, N, 1, fv, fvv);

    sum.resize(M*N);
    for (int i = 0; i < M*N; ++i) {
	sum[i] = (1+fv[i]*fv[i])*fuu[i]
	    + 2*fu[i]*fv[i]*fuv[i]
	    + (1+fu[i]*fu[i])*fvv[i];
	int c = i%M;
	int r = i/M;
	if ((r > 1) && (r < M-2) && (c > 1) && (c < N-2)) {
	    cout << sum[i] << " ";
	}
	if ((i+1)%M == 0) cout << endl;
    }
}
