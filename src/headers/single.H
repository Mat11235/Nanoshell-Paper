/*
 * This file is part of the Quasi-Static Time-Dynamic Project.
 * 
 * Copyright (C) 2025 Alessandro Veltri
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

inline std::complex<double>* pcfc0(std::complex <double> eps_inf, double eps_b, double eps_s=0, double rho=0){
      auto* p = new std::complex<double>[4];
      p[0] = -1./(eps_inf+2.*eps_b);
      p[1] = 2.*p[0];
      p[2] = p[0];
      p[3] = -(eps_inf-eps_b)*p[0];
      return p;
    }

inline std::complex<double>* pcfc1(std::complex <double> eps_inf, double eps_b, double eps_s=0, double rho=0){
      std::complex<double>* p = nullptr;
      return p;
    }

inline std::complex<double>* pcfc2(std::complex <double> eps_inf, double eps_b, double eps_s=0, double rho=0){
      std::complex<double>* p = nullptr;
      return p;
    }
    
inline std::complex<double>* pcfc3(std::complex <double> eps_inf, double eps_b, double eps_s=0, double rho=0){
      std::complex<double>* p = nullptr;
      return p;
    }

    
std::complex<double>** coefficients(std::complex<double> OmeH, std::complex<double> OmeP, std::complex<double> GG, std::complex<double> GamP, std::complex<double>* p0, std::complex<double>* p1 = nullptr, std::complex<double>* p2 = nullptr, int J=0){
      std::complex<double>** A = nullptr;
      A = new std::complex<double>*[3];
      for (int j = 0; j <= 2; j++)  A[j] = new std::complex<double>[3];

      A[0][0] =  OmeH;
      A[0][1] =  0;
      A[0][2] =  0;
      
      A[1][0] =  GG*p0[0];
      A[1][1] =  GG*p0[1] + OmeH;
      A[1][2] =  GG*p0[2];

      A[2][0] =  GamP*p0[0];
      A[2][1] =  GamP*p0[1];
      A[2][2] =  GamP*p0[2] + OmeP;
    
      return A;
    }

std::complex<double>* inhomogeneous(std::complex<double> GG, std::complex<double> GamP, double E0, std::complex<double>* p0, std::complex<double>* p1 = nullptr, std::complex<double>* p2 = nullptr){
      std::complex<double>* B = 0;
      B = new std::complex<double>[3];
      B[0] = GG*E0;
      B[1] = GG*p0[3]*E0;
      B[2] = GamP*(p0[3]-1.)*E0;
      return B;
    }

std::complex<double> F(std::complex<double> eps1, double eps3=0, double rap3=0){
    std::complex<double> eff;
    eff=-0.5*eps1;
    return eff;
    }

double eqn9(double ome, double ome0, std::complex<double> F, double eps_b, double Dome){
    double eq9;
    eq9=Dome*(eps_b-real(F))-2.*(ome - ome0)*imag(F);
    return eq9;
    }    
    
std::complex<double> polarizability(std::complex<double> eps1, std::complex<double> eps2, double eps3=0, double rho=0){
    std::complex<double> alph=(eps1-eps2)/(eps1+2.*eps2);
    return alph;
    }

std::complex<double> numerical_output(double E0, std::complex<double> *q, std::complex<double>* p ){
  std::complex<double> nout;
  nout = p[0]*q[0] + p[1]*q[1] + p[2]*q[2] + p[3]*E0;  
  return nout;
  }
