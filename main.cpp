//Projet Compression d’image par ACP
//Julien Alamelle & Robin Le Meur


#include <iostream>
#include <vector>
#include <boost/gil/gil_all.hpp>
#include <boost/gil/extension/io/jpeg_io.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <utility>

using namespace boost::numeric::ublas;

const double PRES=pow(10,-13);

vector<double> norm(const vector<double> v) {	//norm   : normalise un vecteur
	vector<double> w = v/norm_2(v);				//entrée  : un vecteur de flotant
	return w;}											//sortie : un vecteur de flotant de norme 1

void preVP(const matrix<double> m, vector<double>* v, double* d) {
	bool testArret;
	vector<double> v1((*v).size());
	double n;
	for(unsigned int i=0;i<(*v).size();i++) (*v)(i)=static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
	testArret=true;
	while(testArret){
		(*d)=norm_2((*v));
		(*v)=norm((*v));
		axpy_prod(m,(*v),v1);
		n=norm_2(v1-(*d)*(*v));
		if(n<=PRES) testArret=false;
		else (*v)=v1;
	}
}

int diag_p(const matrix<double> m, matrix<double>* orth, vector<double>* d, double n) {
	if (m.size1()==m.size2()&&m.size1()==orth->size1()&&m.size1()==orth->size2()&&m.size1()==d->size()) {

	matrix<double> m2(m.size1(),m.size1());
	vector<double> v(m.size1());
	double f;
	double a=0, b=0;
	unsigned int i=0;
	m2=m;
	for(;i<m.size1();i++) b+=m(i,i);
	i=0;
	while(a<n*b) {
		std::cout<< i << "\n";
		preVP(m2,&v,&f);
		(*d)(i)=f;
		a+=f;
		i++;
		for (unsigned int j=0;j<m.size1();j++) (*orth)(j,i)=v(j);
		m2 -= (f*outer_prod(v,v));} //déflation
	return i;
	}else std::cout<< "erreur de dimension\n";return 0;}

void diag_n(const matrix<double> m, matrix<double>* orth, vector<double>* d, unsigned int n) {
	if (m.size1()==m.size2()&&m.size1()==orth->size1()&&n==orth->size2()&&n==d->size()&&n>0&&n<=m.size1()) {

	matrix<double> m2(m.size1(),m.size1());
	vector<double> v(m.size1());
	double f;
	m2=m;
	for (unsigned int i=0;i<n;i++) {
		std::cout<< i << "\n";
		preVP(m2,&v,&f);
		(*d)(i)=f;
		for (unsigned int j=0;j<m.size1();j++) (*orth)(j,i)=v(j);
		m2 -= (f*outer_prod(v,v));} //déflation

	}else std::cout<< "erreur de dimension\n";}

void diag(const matrix<double> m, matrix<double>* orth, vector<double>* d) {
	if (m.size1()==m.size2()&&m.size1()==orth->size1()&&m.size1()==orth->size2()&&m.size1()==d->size()) {

	matrix<double> m2(m.size1(),m.size1());
	vector<double> v(m.size1());
	double f;
	m2=m;
	for (unsigned int i=0;i<m.size1();i++) {
		std::cout<< i << "\n";
		preVP(m2,&v,&f);
		(*d)(i)=f;
		for (unsigned int j=0;j<m.size1();j++) (*orth)(j,i)=v(j);
		m2 -= (f*outer_prod(v,v));} //déflation

	}else std::cout<< "erreur de dimension\n";}

void test(int t){
	matrix<double> m(t,t), n(t,t), o(t,t);
	vector<double> v(t);
	vector<vector<double> > tabv(t);

	for(int i=0;i<t;i++) tabv(i)=v;
	for(int i=0;i<t;i++) for(int j=0;j<t;j++) m(i,j)=0;
	for(int i=0;i<t;i++) {m(i,i)=1+rand()%32;m(i,i)=pow(m(i,i),2);std::cout<< m(i,i) << "\n";}
	for(int i=0;i<t;i++) {
		for(int j=0;j<t;j++) (tabv(i))(j)=static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		for(int j=0;j<i;j++) tabv(i) -= tabv(j)*inner_prod(tabv(j),tabv(i));
		tabv(i)=norm(tabv(i));
		for(int j=0;j<t;j++) n(j,i)=(tabv(i))(j);
	}
	//std::cout<< n << "\n";
	axpy_prod(n,m,o);
	axpy_prod(o,trans(n),m);
	diag(m,&n,&v);
	std::cout<< v << "\n";
	//std::cout<< n << "\n";

}

void test_n(int t, unsigned int u){
	matrix<double> m(t,t), n(t,t), o(t,t), p(t,u);
	vector<double> v(u), w(t);
	vector<vector<double> > tabv(t);

	for(int i=0;i<t;i++) tabv(i)=w;
	for(int i=0;i<t;i++) for(int j=0;j<t;j++) m(i,j)=0;
	for(int i=0;i<t;i++) {m(i,i)=1+rand()%32;m(i,i)=pow(m(i,i),2);std::cout<< m(i,i) << "\n";}
	for(int i=0;i<t;i++) {
		for(int j=0;j<t;j++) (tabv(i))(j)=static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		for(int j=0;j<i;j++) tabv(i) -= tabv(j)*inner_prod(tabv(j),tabv(i));
		tabv(i)=norm(tabv(i));
		for(int j=0;j<t;j++) n(j,i)=(tabv(i))(j);
	}
	//std::cout<< n << "\n";
	axpy_prod(n,m,o);
	axpy_prod(o,trans(n),m);
	diag_n(m,&p,&v,u);
	std::cout<< v << "\n";
	//std::cout<< p << "\n";

}

void test_p(int t, double u){
	matrix<double> m(t,t), n(t,t), o(t,t);
	vector<double> v(t);
	vector<vector<double> > tabv(t);
	int x;

	for(int i=0;i<t;i++) tabv(i)=v;
	for(int i=0;i<t;i++) for(int j=0;j<t;j++) m(i,j)=0;
	for(int i=0;i<t;i++) {m(i,i)=1+rand()%32;m(i,i)=pow(m(i,i),2);std::cout<< m(i,i) << "\n";}
	for(int i=0;i<t;i++) {
		for(int j=0;j<t;j++) (tabv(i))(j)=static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
		for(int j=0;j<i;j++) tabv(i) -= tabv(j)*inner_prod(tabv(j),tabv(i));
		tabv(i)=norm(tabv(i));
		for(int j=0;j<t;j++) n(j,i)=(tabv(i))(j);
	}
	//std::cout<< n << "\n";
	axpy_prod(n,m,o);
	axpy_prod(o,trans(n),m);
	x=diag_p(m,&n,&v,u);
	for(int i=0;i<x;i++) std::cout<< v(i) << "\n";
	//std::cout<< n << "\n";

}

pair<matrix<double>, matrix<double> > centrer_reduire(matrix<double> mat, double taux) {
	matrix<double> inverse(mat.size2(),mat.size2()), sym(mat.size1(),mat.size1()), orth(mat.size1(),mat.size1()), matcr(mat.size1(), mat.size2());
	vector<double> moy(mat.size2()), ecart(mat.size2()), temp(mat.size1());
	for (unsigned int i=0;i<mat.size1();i++) temp(i)=1;
	for (unsigned int j=0;j<mat.size2();j++){ moy(j)=0; ecart(j)=0;}
	for (unsigned int i=0;i<mat.size1();i++) {
		for (unsigned int j=0;j<mat.size2();j++) {
			moy(j)+=mat(i,j);
		}
	}
	moy/=mat.size1();
	mat-=outer_prod(temp,moy);

	for (unsigned int i=0;i<mat.size1();i++) {
		for (unsigned int j=0;j<mat.size2();j++) {
			ecart(j)+=mat(i,j)*mat(i,j);
		}
	}
	for (unsigned int i=0;i<mat.size2();i++) for (unsigned int j=0;j<mat.size2();j++) inverse(i,j)=0;
	for (unsigned int j=0;j<mat.size2();j++) inverse(j,j)=1/sqrt(ecart(j));
	axpy_prod(mat,inverse,matcr);
	mat+=outer_prod(temp,moy);
	axpy_prod(matcr,trans(matcr),sym);
	int r = diag_p(sym, orth, temp, taux);
	matrix<double> return1(mat.size1(), r), return2(r, mat.size2());
	for(int i=0; i<temp.size(); i++) temp(i)=(i<r ? 1 : 0);
	axpy_prod(orth, temp, return1);
	axpy_prod(trans(return1), mat, return2);
	return make_pair(return1, return2);

}

template <typename SrcView>
std::vector<matrix<double> > getmat(const SrcView& src) {
    std::vector<matrix<double> > res;
    int nbdecouleur = boost::gil::num_channels<SrcView>::value;
    int nbl = src.height();
    int nbc = src.width();

    for (int i=0; i<nbdecouleur; i++) {
      res.push_back(matrix<double>(nbl, nbc));
    }


    for (int y=0; y<nbl; ++y) {
        typename SrcView::x_iterator src_it = src.row_begin(y);
        for (int x=0; x<nbc; ++x)
            for (int n=0; n<nbdecouleur; ++n)
                res[n](x,y) = src_it[x][n];
    }

    return res;
}

int main()
{
    boost::gil::rgb8_image_t toto;
    double a;
    string s;

    std::vector<matrix<double> > matrices;
    std::pair<matrix<double>, matrix<double> > m;

    /// construction su sourceview toto

    std::cout<< "Choisissez un fichier de format .jpg" << std::endl;
    std::cin>> s;

    boost::gil::jpeg_read_image(s,toto);
    matrices = getmat(boost::gil::const_view(toto));

    std::cout<< "choisissez un taux de conservation" << std::endl;
	std::cin>> a;

    for (int i=0;i<3;i++){
        matrix<double> temp = matrices[i];
        m = centrer_reduire(temp,a);
        axpy_prod(std::get<0>(m), std::get<1>(m), matrices[i]);
    }

    boost::gil::rgb8_image_t titi (matrices[1].size1(), matrices[1].size2());
    boost::gil::rgb8_image_t::view_t v = view(titi);
    for (unsigned int x=0; x<matrices[0].size1(); x++)
    {
      for (unsigned int y=0; y<matrices[0].size2(); y++)
        {
	  v(x,y) = boost::gil::rgb8_pixel_t(matrices[0](x,y),matrices[1](x,y),matrices[2](x,y));
        }
    }

    std::cout<< "Choisissez le nom du fichier, incluez le '.jpg'" << std::endl;
    std::cin>> s;
    boost::gil::jpeg_write_view(s,v);

    return 0;
}
