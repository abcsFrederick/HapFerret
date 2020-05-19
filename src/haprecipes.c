#include <math.h>
#include <stdio.h>
#include "NR.H"
#include "NRUTIL.H"
#include "NRUTIL.C"
#include "hap_hdrs.h"
#include "hap_declare.h"


#include <math.h>
#define NRANSI
#include "nrutil.h"
//#define NMAX 2000
#define NMAX 2000000
#define GET_PSUM \
					for (j=1;j<=ndim;j++) {\
					for (sum=0.0,i=1;i<=mpts;i++) sum += p[i][j];\
					psum[j]=sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}

void amoeba(float **p, float y[], int ndim, float ftol,
	float (*funk)(float []), int *nfunk)
{
	float amotry(float **p, float y[], float psum[], int ndim,
		float (*funk)(float []), int ihi, float fac);
	int i,ihi,ilo,inhi,j,mpts=ndim+1;
	float rtol,sum,swap,ysave,ytry,*psum;

	psum=vector(1,ndim);
	*nfunk=0;
	GET_PSUM
	for (;;) {
		ilo=1;
		ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
		for (i=1;i<=mpts;i++) {
			if (y[i] <= y[ilo]) ilo=i;
			if (y[i] > y[ihi]) {
				inhi=ihi;
				ihi=i;
			} else if (y[i] > y[inhi] && i != ihi) inhi=i;
		}
		rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
		if (rtol < ftol) {
			SWAP(y[1],y[ilo])
			for (i=1;i<=ndim;i++) SWAP(p[1][i],p[ilo][i])
			break;
		}
		fprintf(haplotype_log, "in amoeba,  nfunk = %d\n", nfunk);
		if (*nfunk >= NMAX) nrerror("NMAX exceeded");
		*nfunk += 2;
		ytry=amotry(p,y,psum,ndim,funk,ihi,-1.0);
		if (ytry <= y[ilo])
			ytry=amotry(p,y,psum,ndim,funk,ihi,2.0);
		else if (ytry >= y[inhi]) {
			ysave=y[ihi];
			ytry=amotry(p,y,psum,ndim,funk,ihi,0.5);
			if (ytry >= ysave) {
				for (i=1;i<=mpts;i++) {
					if (i != ilo) {
						for (j=1;j<=ndim;j++)
							p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
						y[i]=(*funk)(psum);
					}
				}
				*nfunk += ndim;
				GET_PSUM
			}
		} else --(*nfunk);
	}
	free_vector(psum,1,ndim);
}
#undef SWAP
#undef GET_PSUM
#undef NMAX
#undef NRANSI


#define NRANSI
#include "nrutil.h"

float amotry(float **p, float y[], float psum[], int ndim,
	float (*funk)(float []), int ihi, float fac)
{
	int j;
	float fac1,fac2,ytry,*ptry;

	ptry=vector(1,ndim);
	fac1=(1.0-fac)/ndim;
	fac2=fac1-fac;
	for (j=1;j<=ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
	ytry=(*funk)(ptry);
	fprintf(haplotype_log, "in amotry 1,  -logP = %f\n", ytry);
	fprintf(hap_test_log, "in amotry 1,  -logP = %f\n", ytry);	
	if (ytry < y[ihi]) {
		y[ihi]=ytry;
		for (j=1;j<=ndim;j++) {
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=ptry[j];
		}
	}
	free_vector(ptry,1,ndim);
	return ytry;
}
#undef NRANSI


#include <math.h>
#define NRANSI
#include "nrutil.h"
#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

float brent(float ax, float bx, float cx, float (*f)(float), float tol,
	float *xmin)
{
	int iter;
	float a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	float e=0.0;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x);
	fprintf(haplotype_log, "in brent frst f call, -logP = %f\n", fw);
	fprintf(hap_test_log, "in brent frst f call, -logP = %f\n", fw);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(u);
		fprintf(haplotype_log, "in brent, its = %d, -logP = %f\n", iter, fu);
		fprintf(hap_test_log, "in brent, its = %d, -logP = %f\n", iter, fu);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	nrerror("Too many iterations in brent");
	*xmin=x;
	return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef NRANSI
#include <math.h>
#define NRANSI
#include "nrutil.h"
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc,
	float (*func)(float))
{
	float ulim,u,r,q,fu,dum;
	int timesthru = 0;

	*fa=(*func)(*ax);
	fprintf(haplotype_log, "in mnbrak 1, *fa func call, x = %f func = %f \n", *ax, *fa);
	fprintf(hap_test_log, "in mnbrak 1, *fa func call, x = %f func = %f \n", *ax, *fa);
	*fb=(*func)(*bx);
	fprintf(haplotype_log, "in mnbrak 1, *fb func call, x = %f func = %f \n", *bx, *fb);
	fprintf(hap_test_log, "in mnbrak 1, *fb func call, x = %f func = %f \n", *bx, *fb);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
	fprintf(haplotype_log, "in mnbrak 1, *fc func call, x = %f func = %f \n", *cx, *fc);
	fprintf(hap_test_log, "in mnbrak 1, *fc func call, x = %f func = %f \n", *cx, *fc);
	while (*fb > *fc) {
		++timesthru;
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u);
			fprintf(haplotype_log, "in mnbrak 1, fu func call, x = %f func = %f \n", u, fu);
			fprintf(hap_test_log, "in mnbrak 1, fu func call, x = %f func = %f \n", u, fu);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
			fprintf(haplotype_log, "in mnbrak 1, fu func call 2, x = %f func = %f \n", u, fu);
			fprintf(hap_test_log, "in mnbrak 1, fu func call 2, x = %f func = %f \n", u, fu);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u);
			fprintf(haplotype_log, "in mnbrak 1, fu func call 3, x = %f func = %f \n", u, fu);
			fprintf(hap_test_log, "in mnbrak 1, fu func call 3, x = %f func = %f \n", u, fu);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
#undef NRANSI

#define NRANSI
#include "nrutil.h"
#define TOL 2.0e-4

int ncom;
float *pcom,*xicom,(*nrfunc)(float []);

void linmin(float p[], float xi[], int n, float *fret, float (*func)(float []))
{
	float brent(float ax, float bx, float cx,
		float (*f)(float), float tol, float *xmin);
	float f1dim(float x);
	void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb,
		float *fc, float (*func)(float));
	int j;
	float xx,xmin,fx,fb,fa,bx,ax;

	ncom=n;
	pcom=vector(1,n);
	xicom=vector(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	// initial guess for brackets:
	ax=0.0;
	xx=0.001;
	//xx=1.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	*fret=brent(ax,xx,bx,f1dim,TOL,&xmin);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free_vector(xicom,1,n);
	free_vector(pcom,1,n);
}
#undef TOL
#undef NRANSI

#include <math.h>
#define NRANSI
#include "nrutil.h"
//#define ITMAX 200
#define ITMAX 1000
#define EPS 1.0e-10
#define FREEALL free_vector(xi,1,n);free_vector(h,1,n);free_vector(g,1,n);

void frprmn(float p[], int n, float ftol, int *iter, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float []))
{
	void linmin(float p[], float xi[], int n, float *fret,
		float (*func)(float []));
	int j,its;
	float gg,gam,fp,dgg;
	float *g,*h,*xi;

	g=vector(1,n);
	h=vector(1,n);
	xi=vector(1,n);
	fp=(*func)(p);
	(*dfunc)(p,xi);
	for (j=1;j<=n;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		fprintf(haplotype_log, "in frprmn 1, its = %d, -logP = %f\n", its, fp);
		fprintf(hap_test_log, "in frprmn 1, its = %d, -logP = %f\n", its, fp);
		linmin(p,xi,n,fret,func);
		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
			FREEALL
			return;
		}
		fp=(*func)(p);
		(*dfunc)(p,xi);
		dgg=gg=0.0;
		for (j=1;j<=n;j++) {
			gg += g[j]*g[j];
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
			FREEALL
			return;
		}
		gam=dgg/gg;
		for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
		fprintf(haplotype_log, "in frprmn 2, its = %d, -logP = %f\n", its, fp);
		fprintf(hap_test_log, "in frprmn 2, its = %d, -logP = %f\n", its, fp);
	}
	nrerror("Too many iterations in frprmn");
}
#undef ITMAX
#undef EPS
#undef FREEALL
#undef NRANSI

#define NRANSI
#include "nrutil.h"

extern int ncom;
extern float *pcom,*xicom,(*nrfunc)(float []);

float f1dim(float x)
{
	int j;
	float f,*xt;

	xt=vector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt);
	free_vector(xt,1,ncom);
	return f;
}
#undef NRANSI


#include <math.h>
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

int Nhaps;
float last_calc_freq, hap_logP, hap_P;

void gcf(float *gammcf, float a, float x, float *gln)
{
	float gammln(float xx);
	void nrerror(char error_text[]);
	int i;
	float an,b,c,d,del,h;

	*gln=gammln(a);
	b=x+1.0-a;
	c=1.0/FPMIN;
	d=1.0/b;
	h=d;
	for (i=1;i<=ITMAX;i++) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=b+an/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
	*gammcf=exp(-x+a*log(x)-(*gln))*h;
}
#undef ITMAX
#undef EPS
#undef FPMIN

#include <math.h>
#define ITMAX 100
#define EPS 3.0e-7

void gser(float *gamser, float a, float x, float *gln)
{
	float gammln(float xx);
	void nrerror(char error_text[]);
	int n;
	float sum,del,ap;

	*gln=gammln(a);
	if (x <= 0.0) {
		if (x < 0.0) nrerror("x less than 0 in routine gser");
		*gamser=0.0;
		return;
	} else {
		ap=a;
		del=sum=1.0/a;
		for (n=1;n<=ITMAX;n++) {
			++ap;
			del *= x/ap;
			sum += del;
			if (fabs(del) < fabs(sum)*EPS) {
				*gamser=sum*exp(-x+a*log(x)-(*gln));
				return;
			}
		}
		nrerror("a too large, ITMAX too small in routine gser");
		return;
	}
}
#undef ITMAX
#undef EPS
float gammq(float a, float x)
{
	void gcf(float *gammcf, float a, float x, float *gln);
	void gser(float *gamser, float a, float x, float *gln);
	void nrerror(char error_text[]);
	float gamser,gammcf,gln;

	if (x < 0.0 || a <= 0.0) nrerror("Invalid arguments in routine gammq");
	if (x < (a+1.0)) {
		gser(&gamser,a,x,&gln);
		return 1.0-gamser;
	} else {
		gcf(&gammcf,a,x,&gln);
		return gammcf;
	}
}

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
	int j;
	long k;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0 || !iy) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ;
			*idum=IA*(*idum-k*IQ)-IR*k;
			if (*idum < 0) *idum += IM;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k;
	if (*idum < 0) *idum += IM;
	j=iy/NDIV;
	iy=iv[j];
	iv[j] = *idum;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

float bico(int n, int k)
{
	float factln(int n);

	return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}

float logbico(int n, int k)
{
	float factln(int n);
	float answer;
	/*printf("logbico n = %d, k = %d \n", n, k);*/
	answer =  factln(n)-factln(k)-factln(n-k);
	return answer;
	
}

float gammln(float xx)
{
	float x,y,tmp,ser;
	static float cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}
float factln(int n)
{
	float gammln(float xx);
	void nrerror(char error_text[]);
	static float a[101];
	float answer;

	if (n < 0)
		 nrerror("Negative factorial in routine factln");
	if (n <= 1) return 0.0;
	if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
	else{
		answer =  gammln(n+1.0);
		return answer;
	}
}

#define NRANSI
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define loopM 7
#define NSTACK 50

void indexx0(unsigned long n, float arr[], unsigned long indx[])
// to index with first value 0, but needs more study...
{
	unsigned long i,indxt,ir=n,itemp,j,k,l=1;
	int jstack=0,*istack;
	float a;
    
	istack=ivector(1,NSTACK);
	//for (j=1;j<=n;j++) indx[j]=j;
	for (j=0;j<n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < loopM) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=1;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir])
			}
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir])
			}
			if (arr[indx[l+1]] > arr[indx[l]]) {
				SWAP(indx[l+1],indx[l])
			}
			i=l+1;
			j=ir;
			indxt=indx[l];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--;   while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j])
			}
			indx[l]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			if (jstack > NSTACK) nrerror("NSTACK too small in indexx0.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_ivector(istack,1,NSTACK);
}
#undef loopM
#undef NSTACK
#undef SWAP
#undef NRANSI

#define NRANSI
#define SWAP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define loopM 7
#define NSTACK 50

void indexx(unsigned long n, float arr[], unsigned long indx[])
{
	unsigned long i,indxt,ir=n,itemp,j,k,l=1;
	int jstack=0,*istack;
	float a;
    
	istack=ivector(1,NSTACK);
	for (j=1;j<=n;j++) indx[j]=j;
	for (;;) {
		if (ir-l < loopM) {
			for (j=l+1;j<=ir;j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1;i>=1;i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir])
			}
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir])
			}
			if (arr[indx[l+1]] > arr[indx[l]]) {
				SWAP(indx[l+1],indx[l])
			}
			i=l+1;
			j=ir;
			indxt=indx[l];
			a=arr[indxt];
			for (;;) {
				do i++; while (arr[indx[i]] < a);
				do j--;   while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j])
			}
			indx[l]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			if (jstack > NSTACK) nrerror("NSTACK too small in indexx.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_ivector(istack,1,NSTACK);
}
#undef loopM
#undef NSTACK
#undef SWAP
#undef NRANSI
float erfcc(float x)
{
	float t,z,ans;

	z=fabs(x);
	t=1.0/(1.0+0.5*z);
	ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
		t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
		t*(-0.82215223+t*0.17087277)))))))));
	return x >= 0.0 ? ans : 2.0-ans;
}
#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

float betacf(float a, float b, float x)
{
	void nrerror(char error_text[]);
	int m,m2;
	float aa,c,d,del,h,qab,qam,qap;

	qab=a+b;
	qap=a+1.0;
	qam=a-1.0;
	c=1.0;
	d=1.0-qab*x/qap;
	if (fabs(d) < FPMIN) d=FPMIN;
	d=1.0/d;
	h=d;
	for (m=1;m<=MAXIT;m++) {
		m2=2*m;
		aa=m*(b-m)*x/((qam+m2)*(a+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		h *= d*c;
		aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
		d=1.0+aa*d;
		if (fabs(d) < FPMIN) d=FPMIN;
		c=1.0+aa/c;
		if (fabs(c) < FPMIN) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if (fabs(del-1.0) < EPS) break;
	}
	if (m > MAXIT) nrerror("a or b too big, or MAXIT too small in betacf");
	return h;
}
#undef MAXIT
#undef EPS
#undef FPMIN

#define NRANSI
#include "nrutil.h"
#define ITMAX 100
#define EPS 3.0e-8

float zbrent(float (*func)(float), float x1, float x2, float tol)
{
	int iter;
	float a=x1,b=x2,c=x2,d,e,min1,min2;
	//float fa=(*func)(a),fb=(*func)(b),fc,p,q,r,s,tol1,xm;
	float fa,fb,fc,p,q,r,s,tol1,xm;
	fa=(*func)(a);
	fb=(*func)(b);
	if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
		nrerror("Root must be bracketed in zbrent");
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += SIGN(tol1,xm);
		fb=(*func)(b);
	}
	nrerror("Maximum number of iterations exceeded in zbrent");
	return 0.0;
}
#undef ITMAX
#undef EPS
#undef NRANSI
float betai(float a, float b, float x)
{
	float betacf(float a, float b, float x);
	float gammln(float xx);
	void nrerror(char error_text[]);
	float bt;

	if (x < 0.0 || x > 1.0) nrerror("Bad x in routine betai");
	if (x == 0.0 || x == 1.0) bt=0.0;
	else
		bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
	if (x < (a+1.0)/(a+b+2.0))
		return bt*betacf(a,b,x)/a;
	else
		return 1.0-bt*betacf(b,a,1.0-x)/b;
}

// following is my own function, not num recipe!
float bicon_p_adj(float p, int N, float ndot)
{
	// adjustments: 1) this is p for occurrence more than ndot times + 1/2 p of occurring ndot times;
	// 2) adjusted so that p = 0.5 if ndot = N*p;
	float n_equil, p_equil, p_uncorr, upperdist, a, b;
	if (p == 0.0){
		if(ndot > 0) return 0;
		else return 0.5; // can this ever happen?  Idea is if p = 0 then all dist is at ndot = 0; but we are splitting P for ndot.
	}
	n_equil = N*p;
	p_equil = 0.5*(betai(n_equil, N-n_equil+1, p) + betai(n_equil+1, N-n_equil+2, p)); 
	p_uncorr = 0.5*(betai(ndot, N-ndot+1, p) + betai(ndot+1, N-ndot+2, p));
	//per notes 6/9/05:
	a = (0.5 - p_equil)/(p_equil*(p_equil - 1));
	b = 1. - a;
	upperdist = a*p_uncorr*p_uncorr + b*p_uncorr;
	if ((ndot < n_equil && upperdist < 0.5) || (ndot > n_equil && upperdist > 0.5)){
		printf("wrong center? in bicon_p_adj!");
	}
	return  upperdist;  // if ndot = n_equil this = .5, right?
}
float binom_conf_logp(float p, int N, float n)
{
	float upperdist, logP, p0;
	if (n == 0){
		p0 = pow(1.0 - p, N);
		upperdist = 1 - p0/2.;
	}
	else{
		upperdist = bicon_p_adj(p, N, (float)n);
	}
	if (upperdist > 0.5) upperdist = MIN(upperdist, 0.99999); // smaller offsets will yield faster action
	if (upperdist < 0.5) upperdist = MAX(upperdist, 0.00001);
	logP = log((1. - upperdist)/upperdist);
	fprintf(hap_test_log, "logP %f\n", logP);
	return logP;
}  																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																	
/*float beta_diff(float n_obs)  // need to consider betai as function of n.  Solving bicon_p_adj for n, 
// given N, P, p, retrieves the n observed, we could replace f with n_obs/N
//{
	//return bicon_p_adj(last_calc_freq, Nhaps, n_obs) - hap_logP;
}*/
float logp_diff(float n_obs)  // need to consider betai as function of n.  Solving bicon_p_adj for n, 
// given N, P, p, retrieves the n observed.
{
	float diff;
	diff = binom_conf_logp(last_calc_freq, Nhaps, n_obs) - hap_logP;
	return diff;
}
/*float simple_diff(float n_obs)  // need to consider betai as function of n.  Solving bicon_p_adj for n, 
// given N, P, p, retrieves the n observed.
{
	return betai(n_obs, Nhaps - n_obs+1, last_calc_freq) - hap_P;
}
float simple_call (float p, int N, float n)  // need to consider betai as function of n.  Solving bicon_p_adj for n, 
{
	float upperdist, logP, p0;
	if (n == 0){
		p0 = pow(1.0 - p, N);
		upperdist = 1 - p0; // not splitting the probability for n itself
	}
	else{
		upperdist =  betai(n, N - n+1, p);
	}
	if (upperdist > 0.5) upperdist = MIN(upperdist, 0.99999); // smaller offsets will yield faster action
	if (upperdist < 0.5) upperdist = MAX(upperdist, 0.00001);	
	logP = log((1. - upperdist)/upperdist);
	return logP;
}*/
// for symmetry and comprehensibility here would be better if "logP" were passed as argument 
float adjust_binom_p(float calc_ct, float logP, int Nchrom)
{
	float newcount, expP;
	// set (local) globals to passed values:
	Nhaps = Nchrom;
	last_calc_freq = calc_ct/(float)Nchrom;
	hap_logP = logP;
	expP = exp(logP);
	hap_P =  1./(expP + 1.);	
	// assuming Psum is actual (dampened) p;
	newcount = zbrent(logp_diff, 0.0, Nchrom, 0.000000001); // beta_diff is now a function of n_obs
	//newcount = zbrent(simple_diff, 0.0001, Nchrom, 0.000000001); // beta_diff is now a function of n_obs
	//p_new = zbrent(&beta_diff, 0, 1, 0.0001);   
	if ((hap_logP > 0 && newcount < calc_ct) || (hap_logP < 0 && newcount > calc_ct)){  
		printf("count smaller, but P > 0.5! \n");
	}
	return newcount;
}
