
#include <iostream>
#include<vector>
#include<algorithm>
#include<utility>
#include<cmath>
#include <limits.h>
using namespace std ;
typedef vector<vector<double>> Matrix  ;
class matrix {
     Matrix v ;
      int n = 0  ;
    Matrix subMatrix (Matrix original , int row ,  int col  ) {
        // return matrix[n-1][n-1]  ;
        Matrix sub ;
        for (int i  = 1 ;  i< row ; ++ i ) {
            vector<double>temp;
            for (int j = 0 ; j < row ; ++ j ) {
                if (j==col )continue ;
                temp.push_back(original[i][j]);
            }
           sub.push_back(temp);
        }
        return sub ;
    }
    Matrix subMatrix (int row  , int col) {
        Matrix sub ;
        for (int i =  0 ; i < n ;++ i ) {
            vector<double>temp;
            for (int j =  0 ;j <  n; ++ j ) {
                 if (j!=col && i != row ) {
                     temp.push_back(v[i][j]);
                 }
            }
           if (temp.empty()==false) {
               sub.push_back(temp);
           }
        }
        return sub;
    }
public:
    int size()const {
        return n ;
    }
    matrix(int s):n(s) {
        Matrix t(n , vector<double>(n)) ;
        v = t ;
    }
    matrix (matrix &m) {
       n = m.n ;
        v = m.v ;
    }
    void transpose() {
        for (int i = 0  ; i< n ; ++ i ) {
            for (int j = i+1; j < n ; ++ j ) {
                swap(v[i][j], v[j][i]) ;
            }
        }
    }
  friend istream &operator>>(istream &read, matrix &m) {
        for (int i = 0 ; i <  m.n; ++ i ) {
             for (int j = 0 ; j < m.n ;++ j ) {
                 read >> m.v[i][j] ;
             }
        }
        return read ;
    }
    friend ostream &operator<< (ostream &out , matrix &m ) {
        for (int i= 0 ;  i<m.n ; ++ i ) {
            for (int j = 0 ; j < m.n ; ++ j ) {
                out << m.v[i][j] << ' '  ;
            }
            out << '\n' ;
        }
        return out ;
    }
    Matrix getmatrix ()const {
        return v ;
    }
    matrix &operator+(matrix &m) {
            for (int i = 0 ; i < n; ++  i ) {
                for (int j = 0 ; j < n ;++ j ) {
                    v[i][j] += m.v[i][j] ;
                }
            }
        return *this ;
    }
    matrix &operator-(matrix &m) {
          for (int i = 0 ; i < n; ++ i ) {
              for (int j = 0 ; j < n; ++ j ) {
                  v[i][j] -= m.v[i][j] ;
              }
          }
        return *this ;
    }
    matrix &operator*(matrix &m ) {
        Matrix t(n , vector<double>(n));
        for (int i = 0 ; i < n; ++ i ) {
            for (int  j = 0 ; j < n; ++ j ) {
                for (int k = 0 ; k < n ;++ k ) {
                    t[i][j] += v[i][k] * m.v[k][j] ;
                }
            }
        }
        v = t ;
        return *this ;
    }
    matrix &operator*(const double &val) {
        for (vector<double>& i : v ) {
            for (double & j : i )j*=val ;
        }
        return *this ;
    }
    matrix &operator=(const matrix &m) {
        for (int i = 0 ; i < n ; ++ i ) {
            for (int j = 0  ;j < n; ++ j ) {
                v[i][j] = m.v[i][j] ;
            }
        }
        return *this ;
    }
    matrix &operator *=(const double &val ) {
        for (vector<double> & i : v) {
            for (double & j : i )j*=val ;
        }
        return *this ;
    }
    matrix &operator *= (matrix &m) {
        Matrix t(n , vector<double>(n));
        for (int i = 0 ; i < n; ++ i ) {
            for (int j = 0 ; j < n ;++ j ) {
                for (int k = 0 ; k < n ;++ k ) {
                    t[i][j] += v[i][k] * m.v[k][j] ;
                }
            }
        }
        v =t ;
        return *this ;
    }
    bool operator==(const matrix &m)const {
        for (int i = 0  ;i < n ;++ i ) {
            for (int j = 0  ; j< n; ++ j ) {
                if (v[i][j] != m.v[i][j])return false ;
            }
        }
        return true ;
    }
    bool operator!=(const matrix &m)const {
        for (int i = 0  ;i < n ;++ i ) {
            for (int j = 0  ; j< n; ++ j ) {
                if (v[i][j] != m.v[i][j])return true ;
            }
        }
        return false ;
    }
    long double  get_determinant (const Matrix &m ) {
        long double  det = 0 ;
        if (m.size()==2) {
            return (m[0][0]* m[1][1]) - (m[0][1] * m[1][0]) ;
        }
        if (m.size()==3) {
            return m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                   + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        }
        for (int i = 0 ; i < m.size() ; ++ i ) {
            det+= m[0][i] * pow(-1 , i) * get_determinant(subMatrix(m , m.size() , i));
        }
        return det ;

    }
    vector<int> getEigenValues() {
         vector<int>eig ;
        for (long long  i = 0 ; i < 1000 ; ++i ) {
            Matrix  temp = v ;
             for (int j = 0 ; j < n; ++ j ) {
                 temp[j][j]-=i ;
             }
           long  double det = get_determinant(temp) ;
            if (det==0) {
                eig.push_back(i);
            }
        }
        return eig ;
    }
    void getInverse() {
        Matrix temp(n , vector<double>(n));
        int k = 0 ;
        for (int i = 0 ; i < n; ++  i) {
            for (int j = 0 ; j < n; ++ j ) {
                temp[i][j] = k & 1 ? -1 : 1 ;
                ++ k ;
            }
        }
        double det=  get_determinant(v) ;
        if (det==0) {
            cout << "there is no inverse : " << endl;

        }
        else {
             for (int i = 0 ; i < n; ++ i ) {
                 for (int j = 0 ; j < n; ++ j) {
                    Matrix sub = subMatrix(i , j );
                     long double det2 = get_determinant(sub);
                     temp[i][j] *= det2;
                 }
             }
            for (int i = 0 ; i < n; ++ i ) {
                for (int j = i+1 ;j< n; ++ j )swap(temp[i][j],temp[j][i]) ;
            }
            for (int i = 0 ; i < n; ++ i ) {
                for (int j = 0 ; j< n; ++ j ) {
                    temp[i][j] *=  (1.0/det) ;
                    cout << temp[i][j] << ' ' ;
                }
                cout << endl;
            }
        }

    }
    bool isSymmetric()const {
        for (int i = 0  ;i < n; ++ i ) {
            for (int j = 0 ; j < n; ++ j ) {
                if (v[i][j]!=v[j][i]) {
                    return false ;
                }
            }
        }
        return true ;
    }
    bool isSkewsymmetric()const {
        for (int i = 0 ;  i < n; ++ i ) {
            for (int j = 0 ; j < n; ++ j ) {
                if (v[i][j]!=-v[j][i]) {
                    return false ;
                }
            }
        }
        return true ;
    }
    bool isUpperTriangular()const {
        for (int i = 1 ; i< n; ++ i ) {
            for (int j = 0 ; j < i ; ++ j ) {
                if (v[i][j]!=0)return false ;
            }
        }
        return true;
    }
    bool isLowerTriangular()const {
        for (int i = 0 ; i < n ;++ i ) {
            for (int j = i +1 ; j< n; ++ j ) {
                if (v[i][j]!=0)return false ;
            }
        }
        return true;
    }
    bool isDiagonalMatrix()const {
        for (int i =  0 ;i < n ;++ i ) {
            for (int j = 0 ;  j< n; ++ j ) {
                if (i==j)continue;
                if (v[i][j]!=0)return false ;
            }

        }
        return true ;
    }
    bool isIdentity()const {
        for (int i = 0 ; i < n; ++ i ) {
            if (v[i][i]!=1) {
                return false ;
            }
            for (int j = 0  ; j < n; ++ j ) {
                if (i==j)continue;
                if (v[i][j]!=0) {
                    return false ;
                }
            }
        }
        return true ;
    }
    bool isScalar()const {
       double k =  v[0][0] ;
        Matrix temp =v ;
        for (int i = 0 ; i < n; ++ i  ) {
            for (int j =  0 ;j < n; ++ j ) {
                temp[i][j]*= 1/k   ;
            }
        }
        for (int i = 0 ; i < n; ++ i ) {
              if (temp[i][i]!=1) {
                  return false ;
              }
            for (int j = 0  ; j < n; ++ j ) {
                if (i==j)continue;
                if (temp[i][j]!=0) {
                    return false ;
                }
            }
        }
        return true ;
    }
    bool isStochastic()const {
        // sum of any row || col = 1
        for (int i = 0 ; i < n; ++ i ) {
             double s = 0  ;
             double s2 = 0 ;
            for (int j = 0  ; j< n; ++ j ) {
                s+=v[i][j] ;
                s2+=v[j][i] ;
            }
            if (s==1 || s2==1)return true ;
        }
        return false ;
    }
};
bool isIdempotent( matrix &m ) {
    // a^k = a
    matrix t = m ;
    m*=m ;
    return t==m ;
}
bool isInvoluntary(matrix &m ) {
    // a ^ 2k = i  or a^2k+1 = a
    matrix t = m ;
    t*=t  ;
    matrix t2= m ;
    t2 = (t2*t2) * m ;
    return t.isIdentity() || t2==m  ;
}
bool isNilPotent(matrix &m ) {
    // a^k = zero matrix
    matrix t = m*m ;
    Matrix v = t.getmatrix() ;
    for (int i = 0 ; i < t.size() ; ++ i ) {
        for (int j = 0 ; j < t.size() ; ++  j) {
            if (v[i][j]!=0) {
                return false ;
            }
        }
    }
    return true ;
}
int main() {



    return 0;
}

/*
6 7 8
3 4 5
3 4 5
*/