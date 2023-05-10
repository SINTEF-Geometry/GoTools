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




#include <vector>
#include <string>
#include <fstream>
#include <assert.h>
#include <time.h>
//#include <ctime>
#include "GoTools/geometry/ObjectHeader.h"
#include "GoTools/geometry/GoTools.h"
//#include "GoTools/geometry/GeomObject.h"
#include "GoTools/geometry/Factory.h"
#include "GoTools/geometry/Utils.h"
#include "GoTools/utils/RegistrationUtils.h"
#include "GoTools/utils/ClosestPointUtils.h"
#include "GoTools/utils/timeutils.h"
#if _OPENMP
#include <omp.h>
#endif


//#define GOTOOLS_LOG

using namespace Go;
using namespace std;

// int show_pts = 10;
int show_pts = 0;


ofstream fileout_status_;


class RegisterPointsStatus
{ // The status starts at 0, ends at 100.

public:

    RegisterPointsStatus(int num_points,
			 const vector<int>& reduce_factors,
			 const vector<double>& tolerances,
			 string& of_status_filename) // If filename has length 0 the status is written to cout.
	: of_status_filename_(of_status_filename),
	  allow_last_iter_clear_(false) // We save some time for calculations after the registration.
    {
	reference_time_ = std::time(0);
	num_points_ = num_points;
	iter_fractions_.resize(reduce_factors.size());
	for (size_t ki = 0; ki < reduce_factors.size(); ++ki)
	{
	    iter_fractions_[ki] = 1.0/reduce_factors[ki];
	}
//        std::cout << "reduce_factors.size(): " << reduce_factors.size() << std::endl;
//	assert(reduce_factors.size() == 2);
	num_iter_.resize(reduce_factors.size());
#if 1 // For the Rhino model with 236 sfs.
	// The num_iter_ is used to estimate the completion percentage.
	// With identity_transformation.
	// num_iter_[0] = 40; // This is not the actual number of iterations but just a qualified guess.
	// num_iter_[1] = 3; // This is not the actual number of iterations but just a qualified guess.
	// With GUI-created transformation.
	assert(reduce_factors.size() == 5);
	// @@sbr We are not takoing into account that the first iterations cost more ...
	// We sort of fix this by reducing the number of iterations for the next levels ...
	// Additionally the time spent is not linear wrt the points, it is less than that.
	num_iter_[0] = 100; // This is not the actual number of iterations but just a qualified guess.
	num_iter_[1] = 100; // This is not the actual number of iterations but just a qualified guess.
	num_iter_[2] = 30; // This is not the actual number of iterations but just a qualified guess.
	num_iter_[3] = 10; // This is not the actual number of iterations but just a qualified guess.
	num_iter_[4] = 2;//3; // Qualified guess. 3 for 1e-03. 20 for 1e-06.
#else // For the original raw model with 12 sfs.
	assert(reduce_factors.size() == 2);
	num_iter_[0] = 200; // This is not the actual number of iterations but just a qualified guess.
	num_iter_[1] = 10; // This is not the actual number of iterations but just a qualified guess.
#endif
	tolerances_ = tolerances;
	num_reg_performed_ = 0;
	curr_iter_level_ = 0;
	curr_compl_perc_ = 1; // We are never at 0 or 100 ...

	updateCompletionStatus();
	// of_status_filename_.clear();
	// fileout_status_ << curr_compl_perc_ << endl;
    }

    void increaseIterationLevel()
    {
	if (num_iter_[curr_iter_level_] > 0)
	{
	    if ((curr_iter_level_ < num_iter_.size() - 1) || (allow_last_iter_clear_))
		num_iter_[curr_iter_level_] = 0;
	}
	if (curr_iter_level_ < num_iter_.size() - 1)
	{
	    curr_iter_level_ += 1;
	}
    }

    // To be called when one step is completed at an iteration level.
    void updatePerformedRegisters(bool done = false)
    {
//	std::cout << "num_iter_[curr_iter_level_]: " << num_iter_[curr_iter_level_] << std::endl;
	num_reg_performed_ += num_points_*iter_fractions_[curr_iter_level_];
	--num_iter_[curr_iter_level_];
	if ((curr_iter_level_ == num_iter_.size() - 1) && (num_iter_[curr_iter_level_] == 0) && (!allow_last_iter_clear_))
	    num_iter_[curr_iter_level_] = 1; // We overrun the clear of iterations on the last level.

	if (!done && (num_iter_[curr_iter_level_] == 0))
	{
	    num_iter_[curr_iter_level_] = 1;
	}

	int new_compl_perc = computeCompletionPercentage();
	if (new_compl_perc > curr_compl_perc_)
	{
	    curr_compl_perc_ = new_compl_perc;
//	    fileout_status_.clear();
	    updateCompletionStatus();
	}

	if (done)
	{
	    clearCurrentIterLevel();
	}
    }

    int currentCompletionPercentage()
    {
	return curr_compl_perc_;
    }

    int nextCompletionPercentage()
    {
	int curr_num_reg = num_points_*iter_fractions_[curr_iter_level_];
	int num_rem_reg = numRemainingRegisters();
	int sum_num_reg = num_reg_performed_ + num_rem_reg;
//	int next_num_rem_reg = num_rem_reg - curr_num_reg;
	int next_num_reg_performed = num_reg_performed_ + curr_num_reg;
	int next_compl_perc = floor(100.0*(double)(next_num_reg_performed)/(double)(sum_num_reg));

	return next_compl_perc;
    }

    // Some time consuming operations should be instructed to update the status file.
    // We then need to supply it with a range of percentages. To avoid decreasing
    // percentages in case we must perform an additional iteration we return a
    // restrictive value.
    // It is also useful to have some space below 100 % to divide among the remaining routines.
    int nextCompletionPercentageWithIncreasedIter()
    {
	int curr_num_reg = num_points_*iter_fractions_[curr_iter_level_];
	int next_num_reg_performed_incr = num_reg_performed_ + curr_num_reg;
	int num_rem_reg = numRemainingRegisters();
	int sum_next_num_reg_incr = next_num_reg_performed_incr + num_rem_reg;
	int next_perc_incr = floor(100.0*(double)(next_num_reg_performed_incr)/(double)(sum_next_num_reg_incr));

#ifndef NDEBUG
	{
	    // Printing to screen
	    int sum_num_reg = num_reg_performed_ + num_rem_reg;
	    cout << "\ncurr_compl_perc_: " << curr_compl_perc_ << endl;
	    cout << "curr_num_reg: " << curr_num_reg << endl;
	    cout << "num_reg_performed_: " << num_reg_performed_ << endl;
	    cout << "num_rem_reg: " << num_rem_reg << endl;
	    cout << "sum_num_reg: " << sum_num_reg << endl;

	    int next_num_rem_reg = num_rem_reg - curr_num_reg;
	    int next_num_reg_performed = num_reg_performed_ + curr_num_reg;
	    int sum_next_num_reg = sum_num_reg;
	    int next_perc = floor(100.0*(double)(next_num_reg_performed)/(double)(sum_next_num_reg));
	    cout << "\nnext_perc: " << next_perc << endl;
	    cout << "next_num_reg_performed: " << next_num_reg_performed << endl;
	    cout << "next_num_rem_reg: " << next_num_rem_reg << endl;
	    cout << "sum_next_num_reg: " << sum_next_num_reg << "\n" << endl;
	    cout << "\nnext_perc_incr: " << next_perc_incr << endl;
	    cout << "next_num_reg_performed_incr: " << next_num_reg_performed_incr << endl;
	    cout << "num_rem_reg: " << num_rem_reg << endl;
	    cout << "sum_next_num_reg_incr: " << sum_next_num_reg_incr << "\n" << endl;
	}
#endif

	return next_perc_incr;
    }

    string statusFilename()
    {
	return of_status_filename_;
    }

    bool lastIterOnCurrLevel()
    {
	return (num_iter_[curr_iter_level_] < 2);
    }

    void printTimeStamps()
    {
	std::ofstream fileout_status_timestamps("tmp/status_timestamps.txt");
	for (size_t ki = 0; ki < timestamps_.size(); ++ki)
	{
	    fileout_status_timestamps << timestamps_[ki].first << " " << timestamps_[ki].second << endl;
	}
    }

    std::time_t referenceTime()
    {
	return reference_time_;
    }

private:

    bool allow_last_iter_clear_; // We may set this to false to save some time for final calculations.
    int num_points_; // The number of 3D points.
    vector<double> iter_fractions_; // The fraction of the points we are testing for the iteration levels.
    vector<int> num_iter_; // Num iterations for each level. Just an empirical guess. May include better estimates using
                           // convergence rate in the future.
    vector<double> tolerances_; // Not used at the moment.
    int curr_iter_level_;
    int num_reg_performed_;
    int curr_compl_perc_;
    string of_status_filename_;
    std::time_t reference_time_; // Seconds since 01-Jan-1970, set from constructor.
    vector<pair<int, int> > timestamps_; // We store time for updates to status_file, % completion and seconds since
                                         // reference.

    // Assuming that the content was altered, writing to file.
    void updateCompletionStatus()
    {
        bool write_to_file = (of_status_filename_.size() > 0);
        if (!write_to_file)
        {
#ifdef GOTOOLS_LOG
            std::string log_level("INFO");
            std::string log_identifier("REGISTRATION_COMPLETION");
            cout << "GOTOOLS LOG: " << log_level << " " << log_identifier << " " << curr_compl_perc_ << endl;
#endif
        }
        else
        {
            std::ofstream fileout_status(of_status_filename_.c_str());
#if 0
            std::time_t current_time = time(0);
            int time_diff = current_time - reference_time_;
            cout << "Adding timestamp to status file (% completion & seconds since start): " <<
                curr_compl_perc_ << " " << time_diff << endl;
            timestamps_.push_back(make_pair(curr_compl_perc_, time_diff));
#endif

            bool use_html_formatting = true;
            if (use_html_formatting)
            {
                const int step = 2;
                const double maxWidth = 800.0;
                const int progress = curr_compl_perc_;
                const int relativeProgress = int((progress/100.0) * maxWidth);
                const string description("Performing the point set registration.");
                const string title = "Registration (step " + std::to_string(step) + " of 2)";
                fileout_status << "<html>\n";
                fileout_status << "<head>\n";
                fileout_status << "<title>" << title << "</title>\n";
                fileout_status << "<link href=\"https://api.eu-cloudflow.eu/portal/twopointo/styles/style.css\" rel=\"stylesheet\" type=\"text/css\">\n";
                fileout_status << "</head>\n";
                fileout_status <<"<body style=\"margin: 20px; padding: 20px;\">\n";
                fileout_status << "<h1>" << title << "</h1>\n";
                fileout_status << "<div style=\"border-radius: 5px; border-color: lightblueblue; border-style:dashed; width: " << maxWidth << "px; height: 80px;padding:0; margin: 0; border-width: 3px;\">\n";
                fileout_status << "<div style=\"position: relative; top: -3px; left: -3px; border-radius: 5px; border-color: lightblue; border-style:solid; width: " << relativeProgress << "px; height: 80px;padding:0; margin: 0; border-width: 3px; background-color: lightblue;\">\n";
                fileout_status << "<h1 style=\"margin-left: 20px;\" >" << progress << "%</h1>\n";
                fileout_status << "</div>\n";
                fileout_status << "</div>\n";
                fileout_status << "<h3>" << description << "</h3>";
                fileout_status << "</body>\n";
                fileout_status << "</html>" << std::endl;
            }
            else
            {
                fileout_status << curr_compl_perc_ << endl;
                fileout_status << "Performing the point set registration." << endl;
            }
        }
    }

    // When iteration has converged, we remove the remaining iter in num_iter_.
    void clearCurrentIterLevel()
    {
	if ((curr_iter_level_ < num_iter_.size() - 1) || (allow_last_iter_clear_))
	    num_iter_[curr_iter_level_] = 0;
    }

    int numRemainingRegisters()
    {
	int num_total_reg = 0;
	for (size_t ki = 0; ki < iter_fractions_.size(); ++ki)
	{
	    num_total_reg += num_points_*iter_fractions_[ki]*num_iter_[ki];
	}

	return num_total_reg;
    }

    int computeCompletionPercentage()
    {
	int num_rem_reg = numRemainingRegisters();
	int curr_completion_perc = floor(100.0*(double)(num_reg_performed_)/(double)(num_reg_performed_ + num_rem_reg));
#if 0
	cout << "num_reg_performed_: " << num_reg_performed_ << ", num_rem_reg: " << num_rem_reg << endl;
	cout << "curr_completion_perc: " << curr_completion_perc << endl;
#endif
	if (curr_completion_perc < 1)
	{
	    curr_completion_perc = 1;
	}
	if (curr_completion_perc > 99)
	{
	    curr_completion_perc = 99;
	}
	return curr_completion_perc;
    }
};

Point centerOfMass(const vector<float>& pts)
{
    const int dim = 3;
    const int num_pts = pts.size()/dim;

    Point pt_sum(0.0, 0.0, 0.0);
    for (size_t ki = 0; ki < pts.size(); ki += dim)
    {
	pt_sum[0] += pts[ki];
	pt_sum[1] += pts[ki+1];
	pt_sum[2] += pts[ki+2];
    }

    Point mass_center = pt_sum/(double)num_pts;

    return mass_center;
}

typedef pair<vector<vector<double> >, Point>  transformation_type;

transformation_type currentTransformation;

transformation_type combine(const transformation_type& first_transf, const transformation_type& second_transf)
{
  vector<vector<double> > rot(3);
  for (int i = 0; i < 3; ++i)
    {
      rot[i].resize(3);
      for (int j = 0; j < 3; ++j)
	{
	  double sum = 0.0;
	  for (int k = 0; k < 3; ++k)
	    sum += second_transf.first[i][k] * first_transf.first[k][j];
	  rot[i][j] = sum;
	}
    }

  Point transl(3);
  for (int i = 0; i < 3; ++i)
    {
      double sum = second_transf.second[i];
      for (int j = 0; j < 3; ++j)
	sum += second_transf.first[i][j] * first_transf.second[j];
      transl[i] = sum;
    }

  return transformation_type(rot, transl);
}


vector<Point> floatsToPoints(const vector<float>& pts_in)
{
  vector<Point> result(pts_in.size() / 3);
  for (int i = 0, idx = 0; i < pts_in.size(); i += 3, ++idx)
    result[idx] = Point(pts_in[i], pts_in[i + 1], pts_in[i + 2]);
  return result;
}


vector<Point> floatsToPoints(const vector<float>& pts_in, const transformation_type& transformation)
{
  vector<Point> result(pts_in.size() / 3);
  vector<vector<double> > rotation = transformation.first;
  Point translation = transformation.second;
  for (int i = 0, idx = 0; i < pts_in.size(); i += 3, ++idx)
    {
      if (idx < show_pts)
	{
	  cout << "  *** Applying tansformation on point " << idx << endl;
	  cout << "    before =";
	  for (int j = 0; j < 3; ++j)
	    cout << " " << pts_in[i+j];
	  cout << endl;
	}
      Point p(3);
      for (int j = 0; j < 3; ++j)
	{
	  double sum = translation[j];
	  for (int k = 0; k < 3; ++k)
	    sum += rotation[j][k] * pts_in[i + k];
	  p[j] = sum;
	}
      result[idx] = p;
      if (idx < show_pts)
	{
	  cout << "    after =";
	  for (int j = 0; j < 3; ++j)
	    cout << " " << p[j];
	  cout << endl;
	}
    }
  return result;
}


double transformationL2(const transformation_type& transformation)
{
  vector<vector<double> > rotation = transformation.first;
  Point translation = transformation.second;
  double sum2 = translation.length2();
  for (int i = 0; i < 3; ++i)
    {
      int next_i = (i + 1) % 3;
      double term = 0.5 * (rotation[i][next_i] - rotation[next_i][i]);
      sum2 += term * term;
    }
  return sum2;
}

// Variant of the Frobenius norm, where we expect the rotation to converge towards
// the identity matrix, and include the translation vector to form a 3x4 matrix.
// @@sbr Our case seems to be dominated by the translation vector, hence no difference from the other version.
double transformationL2_v2(const transformation_type& transformation)
{
  vector<vector<double> > rotation = transformation.first;
  Point translation = transformation.second;
  double sum2 = translation.length2();
  for (int i = 0; i < 3; ++i)
    {
      int next_i = (i + 1) % 3;
      double term0 = 1.0 - rotation[i][i];
      double term1 = rotation[i][next_i];
      double term2 = rotation[next_i][i];
      sum2 += term0*term0 + term1*term1 + term2*term2;
    }
  return sum2;
}


// The avgDist uses the squared L2-norm (which places greater weight on points further away).
double avgDist(const vector<float>& pts1, const vector<float>& pts2, const transformation_type& transformation)
{
  double sum2 = 0.0;
  double sum_L1 = 0.0;
  double sum_L2 = 0.0;
  vector<vector<double> > rotation = transformation.first;
  Point translation = transformation.second;
  int nmb_pts = 0;
  for (int i = 0; i < pts1.size(); i += 3, ++nmb_pts)
    {
      double prev_sum = sum2;
      if (nmb_pts < show_pts)
	{
	  cout << "  *** Dist calculations at " << nmb_pts << endl;
	  cout << "    clp =";
	  for (int j = 0; j < 3; ++j)
	    cout << " " << pts1[i+j];
	  cout << endl;
	  cout << "    pt =";
	  for (int j = 0; j < 3; ++j)
	    cout << " " << pts2[i+j];
	  cout << endl;
	  cout << "    T(pt) =";
	  for (int j = 0; j < 3; ++j)
	    {
	      double sum = translation[j];
	      for (int k = 0; k < 3; ++k)
		sum += rotation[j][k] * pts2[i + k];
	      cout << " " << sum;
	    }
	  cout << endl;
	  cout << "    T(pt) - clp =";
	}
      double sum3 = 0.0;
      for (int j = 0; j < 3; ++j)
	{
	  double sum = translation[j] - pts1[i + j];
	  for (int k = 0; k < 3; ++k)
	    sum += rotation[j][k] * pts2[i + k];
	  sum2 += sum * sum;
	  sum3 += sum * sum;
	  if (nmb_pts < show_pts)
	    cout << " " << sum;
	}
//      std::cout << "sqrt(sum3): " << sqrt(sum3) << std::endl;
      double sqrt_sum3 = sqrt(sum3);
//      std::cout << "DEBUG: sqrt_sum3: " << sqrt_sum3 << std::endl;
      sum_L1 += sqrt(sum3);
      sum_L2 += sum3;
      if (nmb_pts < show_pts)
	cout << endl << "    L2 is " << (sum2 - prev_sum) << endl;
    }
//  std::cout << "DEBUG: Mean deviation: " << sum_L1/(double)nmb_pts << ", mean sum_L2: " << sum_L2/(double)nmb_pts << std::endl;
  return sum2 / (double)nmb_pts;
}


void write_transformation_signed_dists(vector<float>& signed_dists,
//    const vector<float>& pts1, const vector<float>& pts2,
				       const transformation_type& transformation,
				       std::ofstream& fileout)
{
  // ofstream out_str_1("transformed_points.txt");
  // ofstream out_str_2("clp_distances.txt");
  // ofstream out_str_3("transfmat_translvec_clp_distances.txt");

  vector<vector<double> > rotation = transformation.first;
  Point translation = transformation.second;

  for (size_t kj = 0; kj < rotation.size(); ++kj)
  {
      for (size_t ki = 0; ki < rotation[kj].size(); ++ki)
      {
	  fileout << rotation[kj][ki] << " ";
      }
      fileout << endl;
  }
  translation.write(fileout);
  fileout << endl;

  // const int dim = 3;
  int num_pts = signed_dists.size();
  fileout << num_pts << endl;
  for (int i = 0; i < signed_dists.size(); ++i)
    {
      // double sum2 = 0.0;
      // for (int j = 0; j < 3; ++j)
      // 	{
      // 	  double sum = translation[j];
      // 	  for (int k = 0; k < 3; ++k)
      // 	    sum += rotation[j][k] * pts2[i + k];
      // 	  // out_str_1 << sum;
      // 	  // if (j < 2)
      // 	  //   out_str_1 << " ";
      // 	  // else
      // 	  //   out_str_1 << endl;
      // 	  sum -= pts1[i + j];
      // 	  sum2 += sum * sum;
      // 	}
//	fileout << sqrt(sum2) << endl;
	fileout << signed_dists[i] << endl;
    }
}


// We write to the ply file format.
void write_transformed_points_signed_dists(const vector<float>& input_points,
					   const vector<float>& signed_dists,
					   const transformation_type& transformation,
					   std::ofstream& fileout)
{
    const int dim = 3;
    const int num_pts = input_points.size()/dim;
    const int result_size = signed_dists.size()/num_pts;
    //std::cout << "result_size: " << result_size << std::endl;
    if (result_size != 1)
    {
        MESSAGE("Warning: Result includes more than just the signed distances. Sf index and params may be included.");
    }
    fileout << "ply\n";
    fileout << "format ascii 1.0\n";
    fileout << "element vertex " << num_pts << "\n";
    fileout << "property float x\n";
    fileout << "property float y\n";
    fileout << "property float z\n";
    fileout << "property float signed_distance\n";
    fileout << "end_header\n";
    vector<double> all_transf_pts;
    all_transf_pts.reserve(num_pts*dim);
    vector<vector<double> > rotation = transformation.first;
    Point translation = transformation.second;
    for (int i = 0; i < num_pts; ++i)
    {
	vector<double> transf_pt(translation.begin(), translation.end());
	for (int j = 0; j < dim; ++j)
      	{
	    for (int k = 0; k < 3; ++k)
	    {
		transf_pt[j] += rotation[j][k] * input_points[i*dim + k];
	    }
        }
        all_transf_pts.insert(all_transf_pts.end(), transf_pt.begin(), transf_pt.end());
      // 	}
//	fileout << sqrt(sum2) << endl;
        if (result_size == 1)
        {
            fileout << transf_pt[0] << " " << transf_pt[1] << " " << transf_pt[2] << " " << signed_dists[i] << "\n";
        }
        else if (result_size == 4)
        {
            fileout << transf_pt[0] << " " << transf_pt[1] << " " << transf_pt[2] << " " << signed_dists[4*i] << " " <<
               signed_dists[4*i+1] << " " << signed_dists[4*i+2] << " " << signed_dists[4*i+3] << "\n";
        }
    }

#if 0
    Go::PointCloud3D pc_3d(all_transf_pts.begin(), num_pts);
    std::ofstream fileout_pc("tmp/transf_pts.g2");
    pc_3d.writeStandardHeader(fileout_pc);
    pc_3d.write(fileout_pc);
#endif

}

void writeTransformationToFile(const transformation_type& transformation, ofstream& fileout)
{
    // Writing to file the final transformation matrix.
    vector<vector<double> > rotation = transformation.first;
    Point translation = transformation.second;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
            fileout << " " << rotation[i][j];
        fileout << "\n";
    }
    for (int i = 0; i < 3; ++i)
    {
        fileout << " " << translation[i];
    }
    fileout << "\n";
}

void dropTransformation(const transformation_type& transformation, string s)
{
  // cout << "Dropping transformation:" << endl;
  cout << s << endl;
  vector<vector<double> > rotation = transformation.first;
  Point translation = transformation.second;
  for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
	cout << "\t" << rotation[i][j];
      cout << "\t\t" << translation[i] << endl;
    }
  // cout << "Done dropping" << endl;
}

void dropRotation(const transformation_type& transformation, string s, std::ofstream& outfile)
{
  // cout << "Dropping transformation:" << endl;
  cout << s << endl;
  vector<vector<double> > rotation = transformation.first;
  for (int i = 0; i < 3; ++i)
  {
      for (int j = 0; j < 3; ++j)
	  cout << "\t" << rotation[i][j];
      cout << endl;
  }
  // cout << "Done dropping" << endl;
}


void dropTranslation(const transformation_type& transformation, string s)
{
  // cout << "Dropping transformation:" << endl;
  cout << s << endl;
  Point translation = transformation.second;
  for (int i = 0; i < 3; ++i)
  {
      cout << "\t" << translation[i];
  }
  cout << endl;
  // cout << "Done dropping" << endl;
}


void registrationIteration(const vector<float>& pts, const shared_ptr<boxStructuring::BoundingBoxStructure>& structure,
			   double changeL2tol,
			   RegisterPointsStatus& reg_pts_status)//,
//			   int status_max, std::ofstream& status)
{
  int nmb_pts = pts.size() / 3;
  for (int i = 0; i < 10000; ++i)
    {
      // cout << "Running iteration " << i << " for point set of size " << nmb_pts << endl;

      int curr_perc = reg_pts_status.currentCompletionPercentage();
      int next_perc = -1;
      if (reg_pts_status.lastIterOnCurrLevel())
      {
	  next_perc = reg_pts_status.nextCompletionPercentageWithIncreasedIter();
      }
      else
      {
	  next_perc = reg_pts_status.nextCompletionPercentage();
      }

#if 0
      shared_ptr<StatusUpdater> reg_upd;
      if (next_perc - curr_perc > 1)
      {
//	  cout << "We should update the perc more often! curr_perc: " << curr_perc << ", next_perc: " << next_perc << endl;
	  reg_upd = shared_ptr<StatusUpdater>(new StatusUpdater(reg_pts_status.referenceTime()));
	  reg_upd->fileout_status_ = reg_pts_status.statusFilename();
	  reg_upd->curr_perc_ = reg_pts_status.currentCompletionPercentage();
	  reg_upd->curr_perc_local_ = reg_upd->curr_perc_;
	  reg_upd->perc_min_ = reg_upd->curr_perc_;
	  reg_upd->perc_max_ = next_perc;
      }
#endif

//      cout << "Finding closest points." << endl;
      // closestPoints() is the routine which consumes almost all the time (inside this function).
      vector<float> clp = closestPoints(pts, structure, currentTransformation.first, currentTransformation.second);//, reg_upd.get());
//      cout << "Computing avg dist." << endl;
//      std::cout << "DEBUG: Computing avg_dist1." << std::endl;
      double avg_dist1 = avgDist(clp, pts, currentTransformation);
//      std::cout << "avg_dist1: " << avg_dist1 << std::endl;

      Point clp_mass_center = centerOfMass(clp);
      Point pts_mass_center = centerOfMass(pts); // @@sbr201509 This should be computed only once.
      Point transl_mass_center = clp_mass_center - pts_mass_center;
#if 0
      std::cout << "clp_mass_center: (" << clp_mass_center[0] << ", " << clp_mass_center[1] << ", " <<
	  clp_mass_center[2] << ")" << std::endl;
      std::cout << "pts_mass_center: (" << pts_mass_center[0] << ", " << pts_mass_center[1] << ", " <<
	  pts_mass_center[2] << ")" << std::endl;
      std::cout << "transl_mass_center: (" << transl_mass_center[0] << ", " << transl_mass_center[1] << ", " <<
	  transl_mass_center[2] << ")" << std::endl;
#endif

      RegistrationInput regParameters;
      vector<Point> clp_p = floatsToPoints(clp);
      vector<Point> pts_p = floatsToPoints(pts, currentTransformation);
      int max_newton_iterations = regParameters.max_newton_iterations_;
//      cout << "Fine registration." << endl;

#if 0 // @@sbr201509 This does not seem to help us, as least not for current case ...
      std::cout << "Translating with transl_mass_center!" << std::endl;
      // Since this is done on the pts_p which is not used afterwards we only need to update the regResult.translation_ with
      // our transl_mass_center vector.
      for (size_t ki = 0; ki < pts_p.size(); ++ki)
      {
	  pts_p[ki] += transl_mass_center;
      }
#endif

      RegistrationResult regResult = fineRegistration(clp_p, pts_p, false, regParameters);

      transformation_type changeTransformation = transformation_type(regResult.rotation_matrix_, regResult.translation_);

#if 0
      std::cout << "Translating with transl_mass_center!" << std::endl;
      changeTransformation.second += transl_mass_center;
#endif

#if 1
      double changeL2 = transformationL2(changeTransformation);
#else
      std::cout << "DEBUG: Trying alternative transformation norm." << std::endl;
      double changeL2 = transformationL2_v2(changeTransformation);
#endif
      currentTransformation = combine(currentTransformation, changeTransformation);
#if 0
      cout << "Avg dist." << endl;
#endif
//      std::cout << "DEBUG: Computing avg_dist2." << std::endl;
      double avg_dist2 = avgDist(clp, pts, currentTransformation);
//      std::cout << "avg_dist2: " << avg_dist2 << std::endl;

      int newton_iterations = regResult.last_newton_iteration_;
      bool reg_OK = (regResult.result_type_ == RegistrationOK) && (newton_iterations < max_newton_iterations);

#if 0
      // current date/time based on current system
      std::time_t now = time(0);
      // convert now to string form
      char* dt = ctime(&now);
      cout << "The local date and time is: " << dt << endl;

      cout << "Nmb pts: " << nmb_pts << " Iter: " << i;
      cout << " Avg clp dist: " << avg_dist1 << " => " << avg_dist2;
      cout << " Transf L2-chg: " << changeL2;
      cout << " Nmb Nwt it: " << (regResult.last_newton_iteration_ + 1);
      cout << " Rotation:";
      for (int i = 0; i < 3; ++i)
	  for (int j = 0; j < 3; ++j)
	cout << " " << currentTransformation.first[i][j];
      cout << " Transl:";
      for (int j = 0; j < 3; ++j)
	cout << " " << currentTransformation.second[j];
      cout << endl;
#endif
      if (!reg_OK)
	{
	  cout << endl << "******* REGISTRATION FAILED!!! *****" << endl;
	  cout << "  Newton method failing reason: ";
	  switch (regResult.result_type_)
	    {
	    case RegistrationOK:
	      cout << "Did not get close enough in maximum number of iterations (" << max_newton_iterations << ")" << endl;
	      break;
	    case TooFewPoints:
	      cout << "To few input points (must be at least 3, was " << clp.size() << ")" << endl;
	      break;
	    case PointSetSizeDiff:
	      cout << "Different size of input point sets (" << clp.size() << " vs " << pts.size() << ")" << endl;
	      break;
	    case SolveFailed:
	      cout << "Solving linear system failed" << endl;
	      break;
	    default:
	      cout << "Unknown reason (should not happen)" << endl;
	    }
	  exit(1);
	}

      if (changeL2 < changeL2tol)
	{
	  if (nmb_pts > 600000)
	    {
#ifdef GOTOOLS_LOG
	      cout << "Dropping results to files" << endl;
#endif
//	      drop_final(clp, pts, currentTransformation);
	    }
	  reg_pts_status.updatePerformedRegisters(true);
	  break;
	}

//      cout << "Updating performed registers." << endl;
      reg_pts_status.updatePerformedRegisters(false);
    }

}


int main( int argc, char* argv[] )
{
  GoTools::init();

  if (argc != 5 && argc != 6)
    {
      cout << "Usage:  " << argv[0] << " <sf_model.g2> <points.txt> <initial_transf.txt> "
	  "<transf_points_signed_dists.ply> (<final_transf.txt>)" << endl;
      //<completion_status.txt>" << endl;
//	  "<final_transf_signed_dists.txt> <completion_status.txt>" << endl;

      return 1;
    }

#if _OPENMP
  std::cout << "Max available threads = " << omp_get_max_threads() << std::endl;
#endif

  ifstream in_surf(argv[1]);
  ifstream in_pts(argv[2]);
  ifstream in_transf(argv[3]);
  ofstream of_result(argv[4]); // Using the ply format.

  const bool include_sf_and_params = false;

#if 0
//  ofstream of_result(argv[4]); // Line #1-#3: Rotation. #4: Translation. #5: # pts. #6: Signed dist 1st pt. #7: Signed dist 2nd ...
  ofstream of_status(argv[5]); // An integer in the set {0, ..., 100}, an estimated percentage of how much work is done.
  string of_status_filename(argv[5]);
  of_status.clear();
  bool use_html_formatting = true;
  if (use_html_formatting)
  {
      const int step = 1;
      const double maxWidth = 800.0;
      const int progress = 0;//curr_compl_perc;
      const int relativeProgress = int((progress/100.0) * maxWidth);
      std::string description = "Preprocessing the surfaces, initializing.";
      const std::string title = "Registration (step " + std::to_string(step) + " of 2)";
      of_status << "<html>\n";
      of_status << "<head>\n";
      of_status << "<title>" << title << "</title>\n";
      of_status << "<link href=\"https://api.eu-cloudflow.eu/portal/twopointo/styles/style.css\" rel=\"stylesheet\" type=\"text/css\">\n";
      of_status << "</head>\n";
      of_status <<"<body style=\"margin: 20px; padding: 20px;\">\n";
      of_status << "<h1>" << title << "</h1>\n";
      of_status << "<div style=\"border-radius: 5px; border-color: lightblueblue; border-style:dashed; width: " << maxWidth << "px; height: 80px;padding:0; margin: 0; border-width: 3px;\">\n";
      of_status << "<div style=\"position: relative; top: -3px; left: -3px; border-radius: 5px; border-color: lightblue; border-style:solid; width: " << relativeProgress << "px; height: 80px;padding:0; margin: 0; border-width: 3px; background-color: lightblue;\">\n";
      of_status << "<h1 style=\"margin-left: 20px;\" >" << progress << "%</h1>\n";
      of_status << "</div>\n";
      of_status << "</div>\n";
      of_status << "<h3>" << description << "</h3>";
      of_status << "</body>\n";
      of_status << "</html>" << std::endl;
  }
  else
  {
      of_status << 0 << endl;
      of_status << "Preprocessing the surfaces, initializing." << std::endl;
  }
#endif

#ifdef GOTOOLS_LOG
  std::string log_level("INFO");
  std::string log_identifier("REGISTRATION_COMPLETION");
  cout << "GOTOOLS LOG: " << log_level << " " << log_identifier << " " << 0 << endl;
#endif

  ObjectHeader header;
  vector<shared_ptr<GeomObject> > surfaces;
  while (!in_surf.eof())
    {
      header.read(in_surf);
      if (header.classType() == Class_PointCloud)
      {   // The surface file may contain 1D point-cloud (tolerance and unit length prefix + name).
	  const int pc_dim = (header.auxdataSize() == 3) ? header.auxdata(0) : 3;
	  assert((pc_dim == 1) || (pc_dim == 3));
	  if (pc_dim == 1)
	  {
	      Go::PointCloud<1> pc_1d;
	      pc_1d.read(in_surf);
	  }
	  else if (pc_dim == 3)
	  {
	      Go::PointCloud3D pc_3d;
	      pc_3d.read(in_surf);
	  }
	  else
	  {
	      MESSAGE("Did not expect this!");
	  }
      }
      else
      {
	  shared_ptr<GeomObject> obj(Factory::createObject(header.classType()));
	  obj->read(in_surf);
	  surfaces.push_back(obj);
      }
      Utils::eatwhite(in_surf);
    }
  in_surf.close();

  vector<float> pts;

  while (!in_pts.eof())
    {
      for (int j = 0; j < 3; ++j)
	{
	  float f;
	  in_pts >> f;
	  pts.push_back(f);
	}
      Utils::eatwhite(in_pts);
    }
  in_pts.close();

  vector<int> reduce_factors;
  vector<double> tolerances;

  bool use_id_reg = atoi(argv[3]) == 1;

  vector<vector<double> > startRotation;
  Point startTranslation;

  // reduce_factors.push_back(1000);
  // tolerances.push_back(1.0e-6);
  int num_scaled_iter = 0;
  vector<int> level_num_iter;
  vector<int> level_scaled_iter;
  const int dim = 3;
  int num_pts = pts.size()/dim;
  vector<int> perc_each_level;

  // std::cout << "of_status_filename: " << of_status_filename << std::endl;
  // PreprocessPointsStatus prep_pts_status(surfaces.size(), of_status_filename);

  startRotation.resize(3);
  for (int kj = 0; kj < 3; ++kj)
  {
      startRotation[kj].resize(3);
      for (int ki = 0; ki < 3; ++ki)
      {
	  in_transf >> startRotation[kj][ki];
      }
  }
  startTranslation.resize(3);
  startTranslation.read(in_transf);

  currentTransformation = transformation_type(startRotation, startTranslation);
  transformation_type initTransformation = transformation_type(startRotation, startTranslation);

#ifdef GOTOOLS_LOG
  dropTransformation(currentTransformation, "  Input rotation and transformation:");
  cout << "Preprocessing surface model ..." << endl;
#endif

  double ts = getCurrentTime();
  shared_ptr<boxStructuring::BoundingBoxStructure> structure;
  try
  {
      structure = preProcessClosestVectors(surfaces, 200.0);//, &of_status_filename);
  }
  catch (...)
  {
      std::cout << "Failed preprocessing the surfaces!" << std::endl;
      return 1;
  }
  double te = getCurrentTime();
//  std::cout << "DEBUG: Done with the preprocessing, time spent: " << te - ts << std::endl; 

#ifdef GOTOOLS_LOG
  cout << "... done" << endl;
#endif

#if 0
  // We write to screen the size of the largest structure box.
  double voxel_length = structure->voxel_length();
  Point big_vox_low = structure->big_vox_low();
  std::cout << "voxel_length: " << voxel_length << ", big_vox_low: " << big_vox_low << std::endl;
#endif

#if 1
  // CloudFlow: Empirically we experience convergence noise starting
  // at x e-06 L2-change, due to floating point precision input pts I
  // guess.
#if 1
  // L2, i.e. squared, so pure translation (no rotation) w/L2-norm 1e-04 => 1e-02 translation.
  // We are not satisfied with less than 1e-02 translation as this will accumulate, should aim for 1e-03.
  // @@sbr Not sure if less than 1e-06 will converge due to floating point precision.
  const double final_tol = 1.0e-04;
  reduce_factors.push_back(10000);
  tolerances.push_back(0.1*final_tol);//6);//1.0e-09);//4.0e-06);//1.0e-05);//4.0e-06);//10);//7);//5);
  reduce_factors.push_back(1000);
  tolerances.push_back(0.1*final_tol);//6);//1.0e-08);//8);//4.0e-06);//1.0e-05);//4.0e-06);//10);//7);//5);
  reduce_factors.push_back(100);
  tolerances.push_back(0.1*final_tol);//6);//1.0e-07);//8);//4.0e-06);//1.0e-05);//4.0e-06);//10);//7);//5);
  reduce_factors.push_back(10);
  tolerances.push_back(0.1*final_tol);//6);//1.0e-06);//2.0e-05);//1.0e-05);//4.0e-06);//1.0e-05);//4.0e-06);//10);//7);//5);
  reduce_factors.push_back(1);
  tolerances.push_back(final_tol);//5);//1.0e-05);//4.0e-05);//1.0e-04);//4.0e-05);//1.0e-05);//4.0e-06);//4);//3);
#else // August 2015.
  reduce_factors.push_back(100);
  tolerances.push_back(1.0e-5);//1.0e-04);//4.0e-05);//1.0e-05);//4.0e-06);//4);//3);
  reduce_factors.push_back(1);
  tolerances.push_back(1.0e-3);//1.0e-04);//4.0e-05);//1.0e-05);//4.0e-06);//4);//3);
#endif
  // We set the minimum points which will define the model.
  const int min_num_sample_pts = 100;//300;//100;
  for (size_t ki = 0; ki < reduce_factors.size(); ++ki)
  {
      int num_sample_pts = num_pts/reduce_factors[ki];
#ifdef GOTOOLS_LOG
      std::cout << "red_factor: " << reduce_factors[ki] << ", num_sample_pts: " << num_sample_pts << std::endl;
#endif
      if (num_sample_pts < min_num_sample_pts)
      {
	  reduce_factors[ki] *= (float)num_sample_pts/(float)min_num_sample_pts;
	  if (reduce_factors[ki] < 1)
	  {
	      reduce_factors[ki] = 1;
	  }
	  num_sample_pts = num_pts/reduce_factors[ki];
#ifdef GOTOOLS_LOG
	  std::cout << "ki: " << ki << ", red_factor: " << reduce_factors[ki] << ", num_sample_pts: " << num_sample_pts << std::endl;
#endif
      }
  }

  std::string empty_filename(""); // Using an empty filename the status is written to cout.
  RegisterPointsStatus reg_pts_status(num_pts, reduce_factors, tolerances, empty_filename);//, of_status_filename);

  // // Based on empirical data ...
  // // The estimated number of iterations on each level. Assuming linearity (wrt to sample size).
  // level_num_iter.push_back(200);
  // level_num_iter.push_back(10);
  // level_scaled_iter.resize(2);
  // perc_each_level.resize(level_num_iter.size());

  // for (int ki = 0; ki < reduce_factors.size(); ++ki)
  // {
  //     level_scaled_iter[ki] = num_pts*level_num_iter[ki]/reduce_factors[ki];
  //     num_scaled_iter += level_scaled_iter[ki];
  // }
  // for (int ki = 0; ki < perc_each_level.size(); ++ki)
  // {
  //     perc_each_level[ki] = 100*level_scaled_iter[ki]/num_scaled_iter;
  //     cout << "Percentage for level " << ki << ": " << perc_each_level[ki] << endl;
  // }

#else
  cout << "These values are just for testing!" << endl;
  reduce_factors.push_back(10000);
  tolerances.push_back(1.0e-5);
  reduce_factors.push_back(1000);
  tolerances.push_back(1.0e-3);
  std::string empty_filename(""); // Using an empty filename the status is written to cout.
  RegisterPointsStatus reg_pts_status(num_pts, reduce_factors, tolerances, empty_filename);//, of_status_filename);
#endif

#ifndef NDEBUG
  {
//      shared_ptr<StatusUpdater> reg_upd_dummy;
      cout << "Computing closest dist for input matrix." << endl;
      vector<float> signed_dists_init = closestSignedDistances(pts, structure,
							       currentTransformation.first, currentTransformation.second);//,
      //reg_upd_dummy.get());
      // cout << "Writing to file." << endl;
      // cout << "clp.size(): " << clp.size() << ", pts.size(): " << pts.size() << endl;
      std::ofstream fileout_debug("tmp/input_mat_signed_dists.txt");
      write_transformation_signed_dists(signed_dists_init,
					currentTransformation,
					fileout_debug);
      cout << "Done writing input dists to file." << endl;
  }
#endif //NDEBUG

  for (int i = 0; i < reduce_factors.size(); ++i)
  {
      int red_fact = reduce_factors[i];
      if (red_fact == 1)
      {
	  double t0 = getCurrentTime();
//	  std::cout << "DEBUG: Starting the full registration." << std::endl;
	  registrationIteration(pts, structure, tolerances[i], reg_pts_status);
	  double t1 = getCurrentTime();
//	  std::cout << "DEBUG: Done with the full registration, time spent: " << t1 - t0 << std::endl; 
     }
      else
      {
	  // This is a rough registration with a subset of the points.
	  vector<float> few_pts;
	  for (int j = 0, idx = 0; j < pts.size(); j += 3, ++idx)
	      if ((idx % red_fact) == 0)
		  for (int k = 0; k < 3; ++k)
		      few_pts.push_back(pts[j + k]);
	  double t0 = getCurrentTime();
//	  std::cout << "DEBUG: Starting the rough registration." << std::endl;
	  registrationIteration(few_pts, structure, tolerances[i], reg_pts_status);
	  double t1 = getCurrentTime();
//	  std::cout << "DEBUG: Done with the rough registration, time spent: " << t1 - t0 << std::endl;
      }
      reg_pts_status.increaseIterationLevel();
  }

#ifdef GOTOOLS_LOG
  dropTransformation(currentTransformation, "  Final rotation and transformation:");
  cout << "Main diagonal entries diff from 1.0:";
  for (int i = 0; i < 3; ++i)
      cout << " " << (currentTransformation.first[i][i] - 1.0);
  cout << endl;

  cout << "Fetching closest points." << endl; // Used to compute signed dists.
#endif

  if (argc == 6)
  {
      // Writing to file the final transformation matrix.
      ofstream of_result_mat(argv[5]);
      writeTransformationToFile(currentTransformation, of_result_mat);
  }

  vector<float> clp;
#if 1
  int curr_perc = reg_pts_status.currentCompletionPercentage();
  int next_perc = 100;
#if 0
  shared_ptr<StatusUpdater> reg_upd;
//  cout << "curr_perc: " << curr_perc << ", next_perc: " << next_perc << endl;
  if (next_perc - curr_perc > 1)
  {
      //    cout << "We should update the perc more often! curr_perc: " << curr_perc << ", next_perc: " << next_perc << endl;
      reg_upd = shared_ptr<StatusUpdater>(new StatusUpdater(reg_pts_status.referenceTime()));
      reg_upd->fileout_status_ = reg_pts_status.statusFilename();
      reg_upd->curr_perc_ = reg_pts_status.currentCompletionPercentage();
      reg_upd->curr_perc_local_ = reg_upd->curr_perc_;
      reg_upd->perc_min_ = reg_upd->curr_perc_;
      reg_upd->perc_max_ = next_perc;
  }
#endif

#if 0
  clp = closestPoints(pts, structure, currentTransformation.first, currentTransformation.second, reg_upd.get());
#else
  double t0 = getCurrentTime();
//  std::cout << "DEBUG: Calculating the signed distance." << std::endl; 

  vector<float> signed_dists;
  if (include_sf_and_params)
  {
      signed_dists = closestSignedDistanceSfParams(pts, structure,
                                                   currentTransformation.first, currentTransformation.second);//,
  }
  else
  {
      signed_dists = closestSignedDistances(pts, structure,
          currentTransformation.first, currentTransformation.second);//,
  }

  //reg_upd.get());
  double t1 = getCurrentTime();
//  std::cout << "DEBUG: Done calculating the signed distance, time spent: " << t1 - t0 << std::endl; 
#endif
#else
  clp.resize(pts.size());
  std::fill(pts.begin(), pts.end(), 0.0);
  cout << "Creating dummy pts to debug the code!" << endl;
#endif

  cout << "Writing to file." << endl;
  //cout << "clp.size(): " << clp.size() << ", pts.size(): " << pts.size() << endl;
#if 0
  write_transformation_signed_dists(signed_dists,
				    currentTransformation,
				    of_result);
#else

#if 0
  std::cout << "Altering the input transformation, setting it to the identity!" << std::endl;
  transformation_type shuffleTransformation = currentTransformation;
  shuffleTransformation.first[0][0] = 0.0;
  shuffleTransformation.first[0][1] = 0.0;
  shuffleTransformation.first[0][2] = 1.0;
  shuffleTransformation.first[1][0] = 1.0;
  shuffleTransformation.first[1][1] = 0.0;
  shuffleTransformation.first[1][2] = 0.0;
  shuffleTransformation.first[2][0] = 0.0;
  shuffleTransformation.first[2][1] = 1.0;
  shuffleTransformation.first[2][2] = 0.0;
  std::cout << "Using shuffle transformation instead!!!" << std::endl;
#endif

  ts = getCurrentTime();
//  std::cout << "DEBUG: Starting writing the points." << std::endl;

#if 0
  MESSAGE("Resetting to input transformation!");
  currentTransformation = initTransformation;
#endif

  write_transformed_points_signed_dists(pts,
					signed_dists,
					currentTransformation,
					of_result);
  te = getCurrentTime();
//  std::cout << "DEBUG: Done writing the points, time spent: " << te - ts << std::endl;
#endif

#if 0
  cout << "Printing timestamps file! Remove when done debugging!" << endl;
  reg_pts_status.printTimeStamps();
#endif

#if 0
  std::time_t current_time = time(0);
  int time_diff = current_time - reg_pts_status.referenceTime();  
  cout << 100 << " " << time_diff << endl;
#endif

  // clear() does not work, we fetch the file once again.
//  of_status.clear();

#if 0
  std::ofstream of_status_final(of_status_filename);
  if (use_html_formatting)
  {
      const int step = 2;
      const double maxWidth = 800.0;
      const int progress = 100;//curr_compl_perc;
      const int relativeProgress = int((progress/100.0) * maxWidth);
      std::string description = "Done with the registration.";
      const std::string title = "Registration (step " + std::to_string(step) + " of 2)";
      of_status_final << "<html>\n";
      of_status_final << "<head>\n";
      of_status_final << "<title>" << title << "</title>\n";
      of_status_final << "<link href=\"https://api.eu-cloudflow.eu/portal/twopointo/styles/style.css\" rel=\"stylesheet\" type=\"text/css\">\n";
      of_status_final << "</head>\n";
      of_status_final <<"<body style=\"margin: 20px; padding: 20px;\">\n";
      of_status_final << "<h1>" << title << "</h1>\n";
      of_status_final << "<div style=\"border-radius: 5px; border-color: lightblueblue; border-style:dashed; width: " << maxWidth << "px; height: 80px;padding:0; margin: 0; border-width: 3px;\">\n";
      of_status_final << "<div style=\"position: relative; top: -3px; left: -3px; border-radius: 5px; border-color: lightblue; border-style:solid; width: " << relativeProgress << "px; height: 80px;padding:0; margin: 0; border-width: 3px; background-color: lightblue;\">\n";
      of_status_final << "<h1 style=\"margin-left: 20px;\" >" << progress << "%</h1>\n";
      of_status_final << "</div>\n";
      of_status_final << "</div>\n";
      of_status_final << "<h3>" << description << "</h3>";
      of_status_final << "</body>\n";
      of_status_final << "</html>" << std::endl;
  }
  else
  {
      of_status_final << 100 << std::endl;
      of_status_final << "Done with the registration." << std::endl;
  }
#endif
  
#ifdef GOTOOLS_LOG
  cout << "GOTOOLS LOG: " << log_level << " " << log_identifier << " " << 100 << endl;
#endif
}

