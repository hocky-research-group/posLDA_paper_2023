/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2021 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Colvar.h"
#include "ActionRegister.h"
#include "tools/Pbc.h"
#include "tools/File.h"           // Input and output from files 
#include "tools/Matrix.h"         // Linear Algebra operations
#include <sstream>
#include <cmath>
#include "tools/Communicator.h"   // All the MPI related stuffs

namespace PLMD {
namespace colvar {

//+PLUMEDOC;
/*
Describe how to use LDA_PROJ.
*/
//+ENDPLUMEDOC


class lda_proj : public Colvar {

private:
  bool pbc, serial;
  //Communicator comm;                           // Object of communicator type (MPI)
  std::string prec_f_name;      		// precision file name
  std::string ref_f_name;       		// reference file name
  std::string vec_f_name; 			// lda vector file
  IFile in_;             			// create an object of class IFile
  Log out_;
  Matrix<double> ref_str;       	        // coords of reference
  Matrix<double> mobile_str;    		// coords of mobile
  Matrix<double> prec;        			// precision data
  Matrix<double> rotation;
  std::vector<double> lda_vec;                  // lda vector

  //Matrix<double> derv_;
  Matrix<double> derv_numeric;
  Matrix<double> mobile_str_copy;
  void readinputs();            		// reads the input data
  double ld;
  std::vector<AtomNumber> atom_list;            // list of atoms
  const double SMALL = 1.0E-30;
  const double delta = 0.00001;
public:
  static void registerKeywords( Keywords& keys );
  explicit lda_proj(const ActionOptions&);
  double determinant(int n, double *B);
  void kabsch_rot_mat();   		// gives rotation matrix
  double cal_lda_proj();    		// calculates the lda projection
  //void grad_lda(double);        	// calculates the gradient  
  void numeric_lda();        		// calculates the numeric gradient
  // active methods:
  void calculate() override;
};

PLUMED_REGISTER_ACTION(lda_proj, "LDA_PROJ")

// register keywords function
void lda_proj::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("compulsory", "PRECISION", "This file contains precision data for the system." );
  keys.add("compulsory", "REFERENCE", "This file contains coordinates of the reference structure.");
  keys.add("atoms","GROUP","the group of atoms we are calculating the dipole moment for");
  keys.add("compulsory", "VECTOR", "The LDA unit vector");
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
}

// constructor function
lda_proj::lda_proj(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  serial(false),
  prec_f_name(""),
  ref_f_name(""),   
  vec_f_name("")   // Note! no comma here in the last line.
{
  parseFlag("SERIAL",serial);
  parseAtomList("GROUP",atom_list);
  parse("REFERENCE", ref_f_name);
  parse("PRECISION", prec_f_name);
  parse("VECTOR", vec_f_name);
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  checkRead();

  log.printf("  of %u atoms\n",static_cast<unsigned>(atom_list.size()));
  for(unsigned int i=0; i<atom_list.size(); ++i) {
    log.printf("  %d", atom_list[i].serial());
  }

  if(pbc) log.printf("\n using periodic boundary conditions\n");
  else log.printf("\n without periodic boundary conditions\n");
  
  addValueWithDerivatives(); setNotPeriodic();

  requestAtoms(atom_list);

  // call the readinputs() function here
  readinputs();
  
}

// read inputs function
void lda_proj::readinputs()
{
	unsigned N=getNumberOfAtoms();
	// read ref coords
	in_.open(ref_f_name);

	ref_str.resize(N,3); prec.resize(N,N);
	derv_numeric.resize(N,3);
	//log.printf("N = %d", N);

	std::string line_, val_;
	unsigned c_=0;

	while (c_ < N)
	{
		in_.getline(line_);
		std::vector<std::string> items_;
		std::stringstream check_(line_);

		while(std::getline(check_, val_, ' ')){ items_.push_back(val_); }
		for(int i=0; i<3; ++i){ ref_str(c_,i) = std::stold(items_[i]); }
		c_ += 1;
	}
	in_.close();

	//read precision
	in_.open(prec_f_name);
	
	std::string line, val;
	unsigned int c = 0;
	
	while(c < N)
	{
		in_.getline(line);
		
		// vector for storing the objects
		std::vector<std::string> items;

                // stringstream helps to treat a string like an ifstream!
		std::stringstream check(line);

                while (std::getline(check, val, ' '))
                {
                        items.push_back(val);
                }

                for(int i=0; i<N; ++i)
                {
                        prec(c, i) = std::stold(items[i]);
                }

                c += 1;

	}
	in_.close();

	// read lda unit vector
	in_.open(vec_f_name);
	unsigned n_=0;
	std::string l_; 
	while (n_ < N*3){ in_.getline(l_); lda_vec.push_back(std::stod(l_)); n_ += 1; }
	lda_vec.resize(N*3);

	/*
	for(unsigned i =0; i<N*3; ++i){std::cout << lda_vec[i] << " ";}
	std::cout << std::endl; 
	*/

	in_.close();

	// open a output file
	//if(comm.Get_rank()==0) out_.open("OUTPUT");
	//out_.open("OUTPUT");

}



double lda_proj::determinant(int n, double *B)
{
   	
   double A[n][n]={0};
   // make a copy first!
   for(int i=0; i<n; ++i){
	   for(int j=0; j<n; ++j){A[i][j] = B[i * n + j];}
   }
   

   //const double SMALL = 1.0E-30;
   
   //  It calculates determinant of a matrix using partial pivoting.

   double det = 1;

   // Row operations for i = 0, ,,,, n - 2 (n-1 not needed)
   for ( int i = 0; i < n - 1; i++ )
   {
      // Partial pivot: find row r below with largest element in column i
      int r = i;
      double maxA = std::abs( A[i][i] );
      for ( int k = i + 1; k < n; k++ )
      {
         double val = std::abs( A[k][i] );
         if ( val > maxA )
         {
            r = k;
            maxA = val;
         }
      }
      if ( r != i )
      {
         for ( int j = i; j < n; j++ ) std::swap( A[i][j], A[r][j] );
         det = -det;
      }

      // Row operations to make upper-triangular
      double pivot = A[i][i];
      if (std::abs( pivot ) < SMALL ) return 0.0;              // Singular matrix

      for ( int r = i + 1; r < n; r++ )                    // On lower rows
      {
         double multiple = A[r][i] / pivot;                // Multiple of row i to clear element in ith column
         for ( int j = i; j < n; j++ ) A[r][j] -= multiple * A[i][j];
      }
      det *= pivot;                                        // Determinant is product of diagonal
   }

   det *= A[n-1][n-1];

   return det;
}

// kabsch rotation 
//double lda_proj::kabsch_rot_mat() {
void lda_proj::kabsch_rot_mat() {
	
  unsigned N=getNumberOfAtoms();

  Matrix<double> mobile_str_T(3,N);
  Matrix<double> prec_dot_ref_str(N,3);
  Matrix<double> correlation(3,3);


  transpose(mobile_str, mobile_str_T);
  mult(prec, ref_str, prec_dot_ref_str);
  mult(mobile_str_T, prec_dot_ref_str, correlation);


  int rw = correlation.nrows();
  int cl = correlation.ncols();
  int sz = rw*cl;

  // SVD part (taking from plu2/src/tools/Matrix.h: pseudoInvert function)

  std::vector<double> da(sz);
  unsigned k=0;

  // Transfer the matrix to the local array
  for (unsigned i=0; i<cl; ++i) for (unsigned j=0; j<rw; ++j) da[k++]=static_cast<double>( correlation(j,i) ); // note! its [j][i] not [i][j]

  int nsv, info, nrows=rw, ncols=cl;
  if(rw>cl) {nsv=cl;} else {nsv=rw;}

  // Create some containers for stuff from single value decomposition
  std::vector<double> S(nsv);
  std::vector<double> U(nrows*nrows);
  std::vector<double> VT(ncols*ncols);
  std::vector<int> iwork(8*nsv);

  // This optimizes the size of the work array used in lapack singular value decomposition
  int lwork=-1;
  std::vector<double> work(1);
  plumed_lapack_dgesdd( "A", &nrows, &ncols, da.data(), &nrows, S.data(), U.data(), &nrows, VT.data(), &ncols, work.data(), &lwork, iwork.data(), &info );
  //if(info!=0) return info;
  if(info!=0) log.printf("info:", info);

  // Retrieve correct sizes for work and rellocate
  lwork=(int) work[0]; work.resize(lwork);

  // This does the singular value decomposition
  plumed_lapack_dgesdd( "A", &nrows, &ncols, da.data(), &nrows, S.data(), U.data(), &nrows, VT.data(), &ncols, work.data(), &lwork, iwork.data(), &info );
  //if(info!=0) return info;
  if(info!=0) log.printf("info:", info);


  // get U and VT in form of traditional 2D array (U_, VT_)
  double U_[nrows][nrows], VT_[ncols][ncols];

  int  c=0;

  for(unsigned int i=0; i<nrows; ++i){ for(unsigned int j=0; j<nrows; ++j){ U_[j][i] = U[c]; c += 1;} } c = 0;  // note! its [j][i] not [i][j]
  for(unsigned int i=0; i<ncols; ++i){ for(unsigned int j=0; j<ncols; ++j){ VT_[j][i] = VT[c]; c += 1;} } c=0; // note! its [j][i] not [i][j]


  // calculate determinants
  double det_u = determinant(nrows, &U_[0][0]);
  double det_vt = determinant(ncols, &VT_[0][0]);

  // check!
  if (det_u * det_vt < 0.0){ for(int i=0; i<nrows; ++i){U_[i][nrows-1] *= -1;} }


  //Matrix<double> rotation(3,3);
  rotation.resize(3,3);
  Matrix<double> u(3,3), vt(3,3);
  for(int i=0; i<3; ++i){ for(int j=0; j<3; ++j){ u(i,j)=U_[i][j]; vt(i,j)=VT_[i][j]; } }

  // get rotation matrix
  mult(u, vt, rotation);

  //return rotation;
}




// calculates lda projection
double lda_proj::cal_lda_proj() {
  //std::cout << "rank calculate = " << comm.Get_rank() << std::endl;
  
  unsigned N=getNumberOfAtoms();
  
  Matrix<double> rotated_obj(N,3);
  // rotate the object
  mult(mobile_str, rotation, rotated_obj);

  //out_.open("OUTPUT");
  //++++++++++++++++++++++++//
  //matrixOut(out_, rotation);
  //++++++++++++++++++++++++//

  // compute the displacement
  std::vector<double> disp(N*3);
  unsigned c=0;
  for(int i=0; i<N; ++i){ for(int j=0; j<3; ++j) { disp[c] = (rotated_obj(i,j)-ref_str(i,j)); c+=1;} }
  
  /*
  // print
  for(unsigned i =0; i<N*3; ++i){std::cout << disp[i] << " ";}
  std::cout << std::endl;
  */
  
  double ld_val = dotProduct(disp, lda_vec);

  return ld_val;
}

// numeric gradient
void lda_proj::numeric_lda(){
        // This function performs numerical derivative.
	unsigned N=getNumberOfAtoms();
	//derv_numeric.resize(N,3);

	unsigned stride;
	unsigned rank;
	if(serial) {
		// when using components the parallelisation do not work
		stride=1;
		rank=0;
	} else {
		stride=comm.Get_size();
		rank=comm.Get_rank();
	}
//	std::cout << "rank = " << rank << std::endl;
//        log.printf("rank=%d \n", rank);

	for(unsigned i=rank; i<N; i+=stride){
		for (unsigned j=0; j<3; ++j){

			mobile_str(i,j) += delta;
			kabsch_rot_mat();
                        derv_numeric(i,j) = ((cal_lda_proj() - ld)/delta);

			mobile_str(i,j) -= delta;
		}

		/*
		// set atom derivatives
	        Vector vi(derv_numeric(i,0), derv_numeric(i,1), derv_numeric(i,2) );
      	        setAtomsDerivatives(i, vi);
		*/
	}

	if(!serial) {
		if(!derv_numeric.getVector().empty()) comm.Sum(&derv_numeric(0,0), derv_numeric.getVector().size());
	}
    

	/*
        for(int atom=0; atom<N; ++atom){
                for(int j=0; j<3; ++j){

                        Matrix<double> delta_vec(N,3); // setting 0.0 to all elements of a matrix
			delta_vec=0.0;
                        delta_vec(atom,j) += delta;

			// changing here!
			for(int a=0; a<N; ++a){for(int b=0; b<3; ++b){ mobile_str(a,b) = delta_vec(a,b) + mobile_str_copy(a,b);} }
			
			kabsch_rot_mat();
			derv_numeric(atom,j) = cal_lda_proj();
                }
        }

        for(int i=0; i<N; ++i){ for(int j=0; j<3; ++j){ derv_numeric(i,j) -= ld; derv_numeric(i,j) /= delta; }}
        */


        // set derivatives (from rank 0)
	//if (rank == 0){
		for(unsigned i=0; i<N; ++i){
			Vector vi(derv_numeric(i,0), derv_numeric(i,1), derv_numeric(i,2) );
			setAtomsDerivatives(i, vi);
		}
	//}


	//++++++++++++++++++++++++++++//
	//if(rank==0) matrixOut(out_, derv_numeric);
	//matrixOut(out_, derv_numeric);
	
	// clear the matrix (very important step!!)
	derv_numeric *= 0;
}


// calculator
void lda_proj::calculate() {

  if(pbc) makeWhole();
  unsigned N=getNumberOfAtoms();

  mobile_str.resize(N,3);
  
  // load the mobile str
  for(unsigned int i=0; i<N; ++i) {
	  Vector pos=getPosition(i);  // const PLMD::Vector
	  for(unsigned j=0; j<3; ++j){
		  mobile_str(i,j) = pos[j];
	  }
  }

  // translating the structure to center of geometry
  double center_of_geometry[3]= {0.0, 0.0, 0.0};
  
  for(unsigned int i=0; i<N; ++i)
  {
	  center_of_geometry[0] += mobile_str(i,0); center_of_geometry[1] += mobile_str(i,1); center_of_geometry[2] += mobile_str(i,2);
  }

  for(int i=0; i<N; ++i)
  { 
	  for(int j=0; j<3; ++j) { mobile_str(i,j) -= (center_of_geometry[j]/N); } 
  }

  kabsch_rot_mat();
  ld = cal_lda_proj();
  //std::cout << "ld=" << ld;
  
  //grad_lda(ld);
  //if(rank==0) out_.open("OUTPUT");
  //matrixOut(out_, derv_);
  
  /*
  // copy the mobile structure
  mobile_str_copy.resize(N,3);
  mobile_str_copy = mobile_str;
  */
  

  //comm.Bcast(mobile_str,0);
  numeric_lda();
/*  
  if(!serial) {
	  if(!derv_numeric.getVector().empty()) comm.Sum(&derv_numeric(0,0), derv_numeric.getVector().size());
  }
*/

  // set derivatives
  /*
  for(unsigned i=0; i<N; ++i){
	  //Vector vi(derv_(i,0), derv_(i,1), derv_(i,2) );
	  Vector vi(derv_numeric(i,0), derv_numeric(i,1), derv_numeric(i,2) );
	  setAtomsDerivatives(i, vi);
  }
  */

  //if(comm.Get_rank()==0) matrixOut(out_, derv_numeric);

  setBoxDerivativesNoPbc();
  setValue(ld);

  
}

}
}



