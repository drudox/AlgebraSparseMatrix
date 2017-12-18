# include "CCSmatrix.H"

using namespace std;
using namespace mg::numeric::algebra ;


int main(){
  
  CCSmatrix<int> cscm1 = {{11,12,13,14,0,0},{0,22,23,0,0,0},{0,0,33,34,35,36},{0,0,0,44,45,0},
                           {0,0,0,0,0,56},{0,0,0,0,0,66}};
  //CCSmatrix<int> cscm1 = {{11,0,0,0,0,0},{12,22,0,0,0,0},{13,23,33,0,0,0},{14,0,34,44,0,0},{0,0,35,45,0,0},
 //                          {0,0,36,0,56,66}} ;  
  cout << "---------------------------------------------------------------------------------------" << endl; 
  CCSmatrix<int> cscm2("input6.dat");
  cout << "---------------------------------------------------------------------------------------" << endl; 
  //cout << cscm2 ;
  cscm2.print() ;
  cout << "---------------------------------------------------------------------------------------" << endl; 
  CCSmatrix<int> cscm3("input4.dat");
  cout << "---------------------------------------------------------------------------------------" << endl; 
  cout << cscm3 ;
  cout << "---------------------------------------------------------------------------------------" << endl; 
 // CCSmatrix<int> cscm4("input3.dat");
 // cout << "---------------------------------------------------------------------------------------" << endl; 
 // cout << cscm4 ;
 // cout << "---------------------------------------------------------------------------------------" << endl; 
  CCSmatrix<double> cscm5("input7.dat");
  cout << "---------------------------------------------------------------------------------------" << endl; 
  cout << cscm5;
  cout << "---------------------------------------------------------------------------------------" << endl; 
  CCSmatrix<int> cscm4("input3.dat");
  cout << "---------------------------------------------------------------------------------------" << endl; 
  //cscm4(2,0,40);
  cout << cscm4 ;
  cout << "---------------------------------------------------------------------------------------" << endl; 
 // cscm4.printCCS();
  cout << "---------------------------------------------------------------------------------------" << endl; 
  CCSmatrix<double> cscm6  = {{1.01, 0, 3.43, 0}, {   0, 4.07, 0, 0}, {0, 0, 0,  3.09}, {1.0, 2.4, 0,0}} ;  
  cout << "---------------------------------------------------------------------------------------" << endl; 
  cout << cscm5 ;
  cout << "---------------------------------------------------------------------------------------" << endl; 
  CCSmatrix<double> cscm8 = {{   0, 0, 1.21, 0}, {3.31,    0, 0, 0}, {0, 0, 0, -6.11}, {1.0, 2.2, 0,0}} ; 
  cout << "---------------------------------------------------------------------------------------" << endl; 
  cout << cscm6 ;
  cout << "---------------------------------------------------------------------------------------" << endl; 
  std::vector<double> v11{1,0,0,0,0,0,0,0} ;
  CCSmatrix<double> cscm7 = {{11,12,0,14,0,0,0,0},{0,22,23,0,25,0,0,0},{31,0,33,34,0,0,0,0},
                              {0,42,0,0,45,46,0,0},{0,0,0,0,55,0,0,0},{0,0,0,0,65,66,67,0},{0,0,0,0,75,0,77,78},
                              {0,0,0,0,0,0,87,88}} ;
  cout << "---------------------------------------------------------------------------------------" << endl; 
  std::vector<double> v22 = cscm7 * v11 ;  
  cout << "---------------------------------------------------------------------------------------" << endl; 
  for(auto & x : v22)
      cout << x << ' ' ;
  cout << endl;     
  cout << "---------------------------------------------------------------------------------------" << endl; 
  //CCSmatrix<double> cscm9 = 
  //cout << cscm9 ;
  cout << "---------------------------------------------------------------------------------------" << endl; 
  CCSmatrix<double> cscm10 = cscm8 * cscm6 ;
  cscm10.print() ;
  //cscm10.printCCS();
  cout << "---------------------------------------------------------------------------------------" << endl; 
  CCSmatrix<double> cscm11("input8.dat");
  CCSmatrix<double> cscm12("input4.dat");
  cout << "---------------------------------------------------------------------------------------" << endl; 
  CCSmatrix<double> cscm13= cscm12*cscm11 ;  
  cout << cscm13;
  cout << "---------------------------------------------------------------------------------------" << endl; 
  CCSmatrix<double> cscm14 = {{11,0,0,14,0,0,0,0},{0,22,0,0,25,0,0,0},{31,0,0,34,0,0,0,0},
                              {0,42,0,0,0,46,0,0},{0,0,0,0,55,0,0,0},{0,0,0,0,65,0,6,0} };
  cout << "---------------------------------------------------------------------------------------" << endl; 
  CCSmatrix<double> cscm15 = {{0,11,0,0,0,1},{1.2,0,0,6,0,0},{0,0,8,3,0,0},{1,0,0,7,0,0},{3,0,0,0,0,0}
                              ,{11,0,0,8,0,0,},{1,0,0,6,9,0},{0,0,0,1,8,0}} ;
  
  cout << "---------------------------------------------------------------------------------------" << endl; 
  CCSmatrix<double> cscm16 = cscm14*cscm15;
  cout << cscm16;
  cout << "---------------------------------------------------------------------------------------" << endl; 
  CCSmatrix<double> cscm17 = cscm15 *cscm14 ;
  //cout << cscm14 ;
  cout << cscm17 ;
  cout << "---------------------------------------------------------------------------------------" << endl; 
  CCSmatrix<int> cscm18("input10.dat");
  //ofstream f("output.dat", ios::out);
  cout << "---------------------------------------------------------------------------------------" << endl; 
  CCSmatrix<int> ccs5("input17.dat");
  std::vector<int> v1 =  {3,4,0,1,6,8,1,19};    
  std::vector<int> v2 =  ccs5*v1 ;
   
  cout << "---------------------------------------------------------------------------------------" << endl; 
   for(auto& x : v2 )
      cout << x << ' ' ;
   std::cout << endl;   
 //  CCSmatrix<double> ccs6("coo_matrix.mtx"); 
 //  cout << ccs6 ;
     
   CCSmatrix<double> ccs7 = cscm11*cscm12;
   cout << ccs7 ;   
  cout << "---------------------------------------------------------------------------------------" << endl; 
   
   
   CCSmatrix<double> ccs16("input4.dat");
   
   std::vector<double> v3 =  {10,3,0,3,1,2}; 
   std::vector<double> v4 =  ccs16*v3 ;
   
   for (auto & x : v4 )
      cout << x << ' ' ;
   cout <<endl;   

   CCSmatrix<double> ccs186("coo_matrix.mtx");
   cout << ccs186 ;   
  cout << "---------------------------------------------------------------------------------------" << endl; 
 
   CCSmatrix<double> ccs188("mat001.mtx");
   cout << ccs188 ;     return 0;      
}





