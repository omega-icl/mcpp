//=============================================================================
//! Real Double-precision Column Vector Class
class dcovector
{
public:
  ///////////////////////////////////////////////
  /////////////////// objects ///////////////////
  ///////////////////////////////////////////////
  CPPL_INT l; //!< vector size
  CPPL_INT cap; //!< vector capacity
  double* array; //!< 1D array to store vector data
  
  ///////////////////////////////////////////////
  ///////////////// constructors ////////////////
  ///////////////////////////////////////////////
  inline dcovector();
  inline dcovector(const dcovector&);
  inline dcovector(const _dcovector&);
  inline dcovector(const CPPL_INT&, const CPPL_INT=0);
  inline dcovector(const char *);
  inline ~dcovector(); //destructor
  
  ///////////////////////////////////////////////
  ////////////////// functions //////////////////
  ///////////////////////////////////////////////
  //////// cast ////////
  inline _zcovector to_zcovector() const;
  
  //////// io ////////
  inline double& operator()(const CPPL_INT&);
  inline double operator()(const CPPL_INT&) const;
  inline dcovector& set(const CPPL_INT&, const double&); //const;
  inline friend std::ostream& operator<<(std::ostream&, const dcovector&);
  inline void write(const char*) const;
  inline void read(const char*);
  
  //////// calc ////////
  inline friend _drovector t(const dcovector&);
  inline friend double nrm2(const dcovector&);
  inline friend CPPL_INT idamax(const dcovector&);
  inline friend double damax(const dcovector&);
  
  //////// misc ////////
  inline void clear();
  inline dcovector& zero();
  inline void chsign();
  inline void copy(const dcovector&);
  inline void shallow_copy(const _dcovector&);
  inline void alias(const dcovector&);
  inline void unalias();
  inline dcovector& resize(const CPPL_INT&, const CPPL_INT=0);
  inline void stretch(const CPPL_INT&);
  inline friend void swap(dcovector&, dcovector&);
  inline friend _dcovector _(dcovector&);
  
  ///////////////////////////////////////////////
  ///////////// numerical operators /////////////
  ///////////////////////////////////////////////
  //////// = ////////
  inline dcovector& operator=(const  dcovector&);
  inline dcovector& operator=(const _dcovector&);
  
  //////// += ////////
  inline dcovector& operator+=(const  dcovector&);
  inline dcovector& operator+=(const _dcovector&);
  
  //////// -= ////////
  inline dcovector& operator-=(const  dcovector&);
  inline dcovector& operator-=(const _dcovector&);
  
  //////// *= ////////
  inline dcovector& operator*=(const     double&);
  
  //////// /= ////////
  inline dcovector& operator/=(const     double&);
  
  //////// unary ////////
  inline friend const dcovector& operator+(const dcovector&);
  inline friend _dcovector operator-(const dcovector&);
  
  //////// + ////////
  inline friend _dcovector operator+(const  dcovector&, const  dcovector&);
  inline friend _dcovector operator+(const  dcovector&, const _dcovector&);
  
  //////// - ////////
  inline friend _dcovector operator-(const  dcovector&, const  dcovector&);
  inline friend _dcovector operator-(const  dcovector&, const _dcovector&);
  
  //////// * ////////
  inline friend _dgematrix operator*(const  dcovector&, const  drovector&);
  inline friend _dgematrix operator*(const  dcovector&, const _drovector&);
  inline friend _dcovector operator*(const  dcovector&, const     double&);
  
  //////// / ////////
  inline friend _dcovector operator/(const  dcovector&, const     double&);
  
  //////// % ////////
  inline friend     double operator%(const  dcovector&, const  dcovector&);
  inline friend     double operator%(const  dcovector&, const _dcovector&);
  
  //////// double ////////
  inline friend _dcovector operator*(const     double&, const  dcovector&);
  
  //////// hadamard ////////
  inline friend _dcovector hadamard(const dcovector&, const dcovector&);
};
