#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

Rcpp::NumericMatrix transpose_counts(Rcpp::NumericMatrix n)
{
  int temp;
  temp=n(1,2); n(1,2)=n(2,1); n(2,1)=temp;
  temp=n(1,3); n(1,3)=n(3,1); n(3,1)=temp;
  temp=n(1,4); n(1,4)=n(4,1); n(4,1)=temp;
  temp=n(2,3); n(2,3)=n(3,2); n(3,2)=temp;
  temp=n(2,4); n(2,4)=n(4,2); n(4,2)=temp;
  temp=n(3,4); n(3,4)=n(4,3); n(4,3)=temp;
  return(n);
}
