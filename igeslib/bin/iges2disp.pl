#! /usr/bin/perl
#
# Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
# Applied Mathematics, Norway.
#
# Contact information: E-mail: tor.dokken@sintef.no                      
# SINTEF ICT, Department of Applied Mathematics,                         
# P.O. Box 124 Blindern,                                                 
# 0314 Oslo, Norway.                                                     
#
# This file is part of GoTools.
#
# GoTools is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version. 
#
# GoTools is distributed in the hope that it will be useful,        
# but WITHOUT ANY WARRANTY; without even the implied warranty of         
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public
# License along with GoTools. If not, see
# <http://www.gnu.org/licenses/>.
#
# In accordance with Section 7(b) of the GNU Affero General Public
# License, a covered work must retain the producer line in every data
# file that is created or manipulated using GoTools.
#
# Other Usage
# You can be released from the requirements of the license by purchasing
# a commercial license. Buying such a license is mandatory as soon as you
# develop commercial activities involving the GoTools library without
# disclosing the source code of your own applications.
#
# This file may be used in accordance with the terms contained in a
# written agreement between you and SINTEF ICT. 
#

open(FILE, "<@ARGV");

sub convert {
   $dpos = index($_[0],D);
   if ( $dpos < 0 )              # Do not convert.
   {
      print $_[0];
   } else                        # Convert
   {
      $beg = substr($_[0],0,$dpos);
      $ext = substr($_[0],$dpos+1,4);
      $result = $beg . "E" . $ext;
      print $result;
   }
}

print "list\n";
while(<FILE>)  # Read line by line from the file
{
   $first4 = substr($_,0,4);     # Get the first 4 characters of the line.

   if ($first4 eq '128,')        # '128,' means a surface.
   {
      @data = split(/,/, $_);    # Read line into a table using , as separator.
      pop(@data);                # Remove the last element in the data array.

# Read in n1, n2, k1, k2.
      $n1 = $data[1]+1;
      $n2 = $data[2]+1;
      $k1 = $data[3]+1;
      $k2 = $data[4]+1;

      print "surf\n";
      print $n1 . " " . $n2 . " " . $k1 . " " . $k2 . " 3 1\n\n";

# Read in et1, the first knot vector.
      $ones = 0;
      $counter = 0;
      $flag = 1;                 # Ready to read first knot vector.
      for ($i = 10; $i < @data; $i++)
      {
	 $counter++;
	 do convert($data[$i]);  # Printing first or second knot vector.
	 if (($counter == ($n1 + $k1)) && ($flag == 1))
	 {
	    $counter = 0;
	    $flag = 2;           # Ready to read second knot vektor.
	    print "\n\n";
	 } else
	 {
	    print " ";
	 }
      }
      print "\n";
      while ($flag == 1)         # Reading first knot vektor.
      {
	 $_ = <FILE>;
         @data = split(/,/, $_); # Read line into a table using , as separator.
	 pop(@data);             # Remove the last element in the data array.
	 for ($i = 0; $i < @data; $i++)
	 {
	    $counter++;
	    if ($flag == 3)      # Skipping ones.
	    {
	       $ones++;
	       $counter++;
	    } else {             # Printing first or second knot vector.
	       do convert( $data[$i] );
            }
	    if (($counter == ($n1 + $k1)) && ($flag == 1))
	    {
	       $counter = 0;
	       $flag = 2;        # Ready to read second knot vektor.
	       print "\n\n";
            } else
            {
	       if (($counter == ($n2 + $k2)) && ($flag == 2))
	       {
	          $counter = 0;
	          $flag = 3;     # Ready to skip ones.
	          print "\n\n";
               } else
               {
	          print " ";
	       }
            }
         }
	 print "\n";
      }

# Read in et2, the second knot vector.
      while ($flag == 2)         # Reading second knot vektor.
      {
	 $_ = <FILE>;
         @data = split(/,/, $_); # Read line into a table using , as separator.
	 pop(@data);             # Remove the last element in the data array.
	 for ($i = 0; $i < @data; $i++)
	 {
	    $counter++;
	    if ($flag == 2)      # Reading second knot vektor.
	    {
	       do convert( $data[$i] );
            } else {             # Skipping ones.
	       $ones++;
	    }
	    if (($counter >= ($n2 + $k2)) && ($flag == 2))
	    {
	       $counter = 0;
	       $flag = 3;        # Ready to skip ones.
	       print "\n\n";
            } else
            {
	       if ($flag == 2)
               {
	          print " ";
	       }
	    }
         }
	 print "\n";
      }

# Remove all the ($n1*$n2) ones.
      while ($flag == 3)
      {
	 $_ = <FILE>;
         @data = split(/,/, $_); # Read line into a table using , as separator.
	 pop(@data);             # Remove the last element in the data array.
	 for ($i = 0; $i < @data; $i++)
	 {
	    $ones++;             # Skipping ones.
	    if ($flag == 4)      # Reading control vertices.
	    {
	       do convert( $data[$i] );
	       $counter++;
	       if (($counter%3) == 0)
	       {
		  print "\n";
	       } else {
		  print " ";
               }
	    }
	    if (($ones == ($n1*$n2)) && ($flag == 3))
	    {
	       $counter = 0;
	       $flag = 4;        # Ready to read control vertices.
	    }
         }
      }


# Read the CV's, one point (x,y,z) for each line.
      while($counter < $n1*$n2*3) # @OBS: should use $flag == 4
      {
	 $_ = <FILE>;
         @data = split(/,/, $_); # Read line into a table using , as seperator.
	 pop(@data);             # Remove the last element in the data array.
	 foreach $key (@data)
	 {
	    if ($counter < $n1*$n2*3)
	    {
 	       do convert( $key );
	    }
	    $counter += 1;
	    if ( ($counter%3) == 0)
	    {
	       print "\n";
            } else {
               print " ";
            }
	 }
      }

      print "\n\n";
   }
}
print "end\n\n";

close(FILE);
