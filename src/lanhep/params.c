#include <math.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <dlfcn.h>

#include "lanhep.h"


typedef double (*extfunc)(double,...);

int P_D_WIDTH=85;
int longest_pdline=0;
		
static List params=0;
static Term rm_zero(Term);
static int legal_expr(Term t);
static List ExtFuncList=0;

extern int opTriHeu;

List all_param_list(void)
	{
	return params;
	}
	
static double cu_rt(double r, double i)
{
	double mod, ph;
/*	printf("cu_rt: %g %g; ",r,i);*/
	mod=sqrt(r*r+i*i);
	ph=atan2(i,r);
/*	printf("mod=%g, ph=%f\n",mod,ph);*/
	mod=pow(mod,1.0/3.0);
	ph/=3;
	return mod*cos(ph);
}

static void proc_param(Term t)
	{
	Atom name, value, comment;
	List li;
	
	if(is_compound(t) && CompoundName(t)==OPR_COMMA)
		{
		Term a,b;
		a=ConsumeCompoundArg(t,1);
		b=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		proc_param(a);
		proc_param(b);
		return;
		}
		
	if(is_atom(t) /*CompoundName(t)!=OPR_EQSIGN*/)
		{ name=t; value=0; comment=0; }
	else 
		{
		Term t1;
		name=ConsumeCompoundArg(t,1);
		t1=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		t=t1;
		if(is_compound(t) && CompoundName(t)==OPR_COLON)
			{
			value=ConsumeCompoundArg(t,1);
			comment=ConsumeCompoundArg(t,2);
			FreeAtomic(t);
			}
		else
			{ value=t; comment=0; }
		}
		

	
	if(!is_atom(name))
		{
		ErrorInfo(101);
		printf(" \'");
		WriteTerm(name);
		printf("\' is not appropriate name for parameter\n");
		FreeAtomic(name);
		FreeAtomic(value);
		FreeAtomic(comment);
		return;
		}
		
	if(comment!=0 && !is_atom(comment))
		{
		ErrorInfo(102);
		printf(" \'");
		WriteTerm(comment);
		printf("\' is not appropriate comment for parameter \'%s\'\n",
				AtomValue(name));
		FreeAtomic(name);
		FreeAtomic(value);
		FreeAtomic(comment);
		return;
		}
		
	value=ProcessAlias(value);
	/*WriteTerm(value);printf(" -> ");*/
	value=rm_zero(value);
	/*WriteTerm(value);puts("");*/
	
	
	if(value !=0 && !legal_expr(value))
		{
		ErrorInfo(103);
		printf(" \'");
		WriteTerm(value);
		printf("\' is not appropriate value for parameter \'%s\'\n",
				AtomValue(name));
		FreeAtomic(name);
		FreeAtomic(value);
		FreeAtomic(comment);
		return;
		}
		
	li=params;
	while(!is_empty_list(li))
		{
		Term t1;
		t1=ListFirst(li);
		if(FunctorName(CompoundFunctor(t1))==name)
			{
			if(value!=0)
				SetCompoundArg(t1,1,value);
			if(comment!=0)
				SetCompoundArg(t1,2,comment);
			
			if(strcmp(AtomValue(name),"Maux"))
				{
				WarningInfo(104);
				printf("Changing parameter \'%s\':",AtomValue(name));
				printf("  value= ");	 WriteTerm(CompoundArg1(t1));
				printf("  comment= ");	 WriteTerm(CompoundArg2(t1));
				puts("");
				}
			return;
			}
		li=ListTail(li);
		}
	
/*	fu=NewFunctor(name,3);*/
	t=MakeCompound3(name,value,comment,0);
/*	SetCompoundArg(t,1,value);
	SetCompoundArg(t,2,comment);*/
	params=AppendLast(params,t);
	ReportRedefined(name,"parameter");
	SetAtomProperty(name,PROP_TYPE,OPR_PARAMETER);
	
	if(opTriHeu)
		tri_reg_prm(name,value);
	
/*	if(IsTermInput())
		{
		printf("Parameter \'%s\' ",AtomValue(name));
		if(value!=0)
			{ printf("   value= ");		WriteTerm(value);}
		if(comment!=0)
			{ printf("   comment= ");	WriteTerm(comment);}
		puts("");
		}*/
	}


Term ProcessParameter(Term t, Term ind)
	{
	Term arg;
	char *s;
	arg=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	if(is_atom(arg))
		{
		s=AtomValue(arg);
		if(s[0]=='?' && s[1]==0)
			{
			List li;
			li=params;
			while(!is_empty_list(li))
				{
				Term t1,t2;
				t1=ListFirst(li);
				printf(AtomValue(FunctorName(CompoundFunctor(t1))));
				t2=CompoundArg1(t1);
				if(t2!=0)
					{ printf(" = "); WriteTerm(t2); }
				t2=CompoundArg2(t1);
				if(t2!=0)
					{ printf(" : "); WriteTerm(t2); }
				if(CompoundArgN(t1,3))
					{ printf(" ( = %f )",FloatValue(CompoundArgN(t1,3))); }
				puts("");
				li=ListTail(li);
				}
			return 0;
			}
		}
	proc_param(arg);
	return 0;
	}

int is_parameter(Atom p)
	{
	if(p==A_SQRT2)
		return 1;
	if(GetAtomProperty(p,A_INFINITESIMAL))
		return 1;
    return (GetAtomProperty(p,PROP_TYPE)==OPR_PARAMETER);
    /*
	li=params;
	while(!is_empty_list(li))
		{
		if(FunctorName(CompoundFunctor(ListFirst(li)))==p)
			return 1;
		li=ListTail(li);
		}
	return 0;
    */
	}


static int legal_expr(Term t)
	{
	char c;
	Functor fu;
	int arity,i,ret;

	ret=1;
	c=AtomicType(t);
	switch(c)
		{
		case 'i':
		case 'f':
			return 1;
		case 'a':
			if(is_parameter(t))
				return 1;
			printf("Error: symbol \'%s\' was not declared.\n",AtomValue(t));
			return 0;
		case 'c':

			fu=CompoundFunctor(t);
			if( (CompoundArity(t)==2 && CompoundName(t)==OPR_PLUS)  ||
				(CompoundArity(t)==1 && CompoundName(t)==OPR_PLUS)  ||
				(CompoundArity(t)==2 && CompoundName(t)==OPR_MINUS) ||
				(CompoundArity(t)==1 && CompoundName(t)==OPR_MINUS) ||
				(CompoundArity(t)==2 && CompoundName(t)==OPR_DIV)   ||
				(CompoundArity(t)==2 && CompoundName(t)==OPR_MLT)   ||
				(CompoundArity(t)==2 && CompoundName(t)==OPR_POW)   ||
				(CompoundArity(t)==1 && CompoundName(t)==A_SIN)     ||
				(CompoundArity(t)==1 && CompoundName(t)==A_COS)     ||
				(CompoundArity(t)==1 && CompoundName(t)==A_ASIN)    ||
				(CompoundArity(t)==1 && CompoundName(t)==A_ACOS)    ||
				(CompoundArity(t)==1 && CompoundName(t)==A_TAN)     ||
				(CompoundArity(t)==1 && CompoundName(t)==A_ATAN)    ||
				(CompoundArity(t)==1 && CompoundName(t)==A_EXP)     ||
				(CompoundArity(t)==2 && CompoundName(t)==A_ATAN2)   ||
				(CompoundArity(t)==3 && CompoundName(t)==A_IF)      ||
				(CompoundArity(t)==1 && CompoundName(t)==A_LOG)     ||
				(CompoundArity(t)==1 && CompoundName(t)==A_SQRT)    ||
				(CompoundArity(t)==1 && CompoundName(t)==A_FABS)    ||
				(GetAtomProperty(CompoundName(t),A_EXT_FUNC)))
			{
			Integer efa=GetAtomProperty(CompoundName(t),A_EXT_FUNC);
			arity=FunctorArity(fu);
			
			if(efa && IntegerValue(efa)!=0  && IntegerValue(efa)!=arity)
				{
				printf("Wrong arguments number in function %s(...).\n",
					AtomValue(CompoundName(t)));
				return 0;
				}
			
			for(i=1;i<=arity;i++)
				{
				if(!legal_expr(CompoundArgN(t,i)))
					ret=0;
				}
			return ret;
			}
			else
			{
				printf("Error: function ");WriteTerm(t);puts(" is unknown.");
				return 0;
			}
		default:
			printf("Internal error: bad expression.\n");
			return 0;
		}
	return 0;
	}


static char wpbuf[128];

static void texWriteParams(int fno, char *name)
	{
	FILE *f1;
	List li;
	if(OutputDirectory!=NULL)
		sprintf(wpbuf,"%s/vars%d.tex",OutputDirectory,fno);
	else
		sprintf(wpbuf,"vars%d.tex",fno);
	f1=fopen(wpbuf,"w");
	if(f1==NULL)
		{
		printf("Can not open file \'%s\' for writing.\n",wpbuf);
		perror("");
		return;
		}
	
	fprintf(f1,"\\begin{tabular}{|l|l|l|} \\hline\n");	
	fprintf(f1,"Parameter & Value & Comment \\\\ \\hline\n");
	li=params;
	while(!is_empty_list(li))
		{
		Term ttt, val, comm;
		ttt=ListFirst(li);
		val=CompoundArg1(ttt);
		comm=CompoundArg2(ttt);
		if(val==0) val=NewInteger(0);
		if(is_integer(val) || is_float(val))
			{
			int sp;
			sp=fprintf(f1,AtomValue(FunctorName(CompoundFunctor(ttt))));
			WriteBlank(f1,6-sp);
			fprintf(f1,"&");
			sp=fWriteTerm(f1,val);
			WriteBlank(f1,20-sp);
			fprintf(f1,"&");
			if(comm!=0)
				{
				if(ListTail(li))
					fprintf(f1,"%s \\\\\n",AtomValue(comm));
				else
					fprintf(f1,"%s ",AtomValue(comm));
				}
			else
				{
				if(ListTail(li))
					fprintf(f1," \\\\\n");
				else
					fprintf(f1," ");
				}
			}
		else
			{
			int sp;
			sp=fprintf(f1,AtomValue(FunctorName(CompoundFunctor(ttt))));
			WriteBlank(f1,6-sp);
			fprintf(f1,"&");
			sp=fWriteTerm(f1,val);
			WriteBlank(f1,20-sp);
			fprintf(f1,"&");
			if(comm!=0)
				{
				if(ListTail(li))
					fprintf(f1,"%s \\\\\n",AtomValue(comm));
				else
					fprintf(f1,"%s ",AtomValue(comm));
				}
			else
				{
				if(ListTail(li))
					fprintf(f1," \\\\\n");
				else
					fprintf(f1," ");
				}
			}
		li=ListTail(li);
		}
	fprintf(f1,"\\\\ \\hline\n\\end{tabular}\n");
	fclose(f1);
	}

extern int EvalPrm;

static void rpl_pow(Term e)
{
	int a,i;
	if(!is_compound(e))
		return;
	a=CompoundArity(e);
	for(i=1;i<=a;i++)
		rpl_pow(CompoundArgN(e,i));
	if(CompoundName(e)==OPR_POW)
		SetCompoundName(e,OPR_CARET);
}

static int has_Q(Term v)
{
	int i;
	int r=0;
	if(is_atom(v) && AtomValue(v)[0]=='Q' && AtomValue(v)[1]==0)
		return 1;
	if(!is_compound(v))
		return 0;
	for(i=1;i<=CompoundArity(v);i++)
		r|=has_Q(CompoundArgN(v,i));
	return r;
}

static int micro_allpr(Atom p, Term val)
{
	char *nm=AtomValue(p);
	if(nm[0]=='B'&&nm[1]=='0')
		return 0;
	if(is_integer(val)||is_float(val))
		return 1;
	return !has_Q(val);
}

List cls_real_matr(void);

void FADeclRealParam(FILE *f)
{
	List li,li2;
	int i;
	fprintf(f,"Scan[ (RealQ[#] = True)&,\n  { "); 
	for(li=params,i=1;li;li=ListTail(li),i++)
	{
		fprintf(f,"%s%c ",AtomValue(CompoundName(ListFirst(li))),
				ListTail(li)?',':' ');
		if(i%10==0) fprintf(f,"\n    ");
	}
	fprintf(f,"} ]\n\n");
	li2=cls_real_matr();
	if(li2==0)
		return;
	
	fprintf(f,"Scan[ (RealQ[#] = True)&,\n  { "); 
	for(li=li2,i=1;li;li=ListTail(li),i++)
	{
		fprintf(f,"%s%c ",AtomValue(ListFirst(li)),
				ListTail(li)?',':' ');
		if(i%8==0) fprintf(f,"\n    ");
	}
	fprintf(f,"} ]\n\n");
	FreeAtomic(li2);
}

void cls_write_decl(FILE *);
void cls_write_matr(FILE *);
extern List fainclude;

void FAWriteParameters(int fno)
{
	char cbuf[128];
	time_t tm;
	List li;
	FILE *f;
	int p;
	int has_GG=0;
	int first=1;
	
	if(OutputDirectory!=NULL)
		sprintf(cbuf,"%s/model%d.h",OutputDirectory,fno);
	else
		sprintf(cbuf,"model%d.h",fno);
	f=fopen(cbuf,"w");
	
	
	fprintf(f,"*     LanHEP output produced at ");
	time(&tm);
	fprintf(f,ctime(&tm));
	if(ModelName)
		fprintf(f,"*     Model named '%s'\n",ModelName);
	fprintf(f,"\n");
	
	fprintf(f,"      double precision Sqrt2, pi, degree, hbar_c2,bogus\n");
	fprintf(f,"      parameter (Sqrt2=1.41421356237309504880168872421D0)\n");
	fprintf(f,"      parameter (pi = 3.1415926535897932384626433832795029D0)\n");
	fprintf(f,"      parameter (degree = pi/180D0)\n");
	fprintf(f,"      parameter (hbar_c2 = 3.8937966D8)\n");
	fprintf(f,"      parameter (bogus = -1D123)\n");
	fprintf(f,"      double complex cI\n      parameter (cI = (0D0, 1D0))\n\n");
	
	fprintf(f,"      double precision Divergence\n      common /renorm/ Divergence\n\n");

		
	p=fprintf(f,"      double precision ");
	
	for(li=params;li;li=ListTail(li))
	{
		Term ttt=CompoundName(ListFirst(li));
		if(strcmp(AtomValue(ttt),"pi")==0)
		  continue;
		if(ttt==A_GG) has_GG=1;
		if(p>60)
		{
			fprintf(f,"\n");
			p=fprintf(f,"      double precision ");
			first=1;
		}
		if(!first)
			p+=fprintf(f,", ");
		else
			first=0;
		p+=fprintf(f,"%s",AtomValue(ttt));
	}
	if(!has_GG)
	{
		if(!first) p+=fprintf(f,", ");
		fprintf(f,"GG");
	}
	fprintf(f,"\n\n");
	
	cls_write_decl(f);
	
	first=1;
	fprintf(f,"      common /mdl_para/\n     &    ");
	p=10;
	for(li=params;li;li=ListTail(li))
	{
		Term ttt=CompoundName(ListFirst(li));
		if(strcmp(AtomValue(ttt),"pi")==0)
		  continue;
		if(p>60)
		{
			fprintf(f,"\n");
			p=fprintf(f,"     &    ");
		}
		p+=fprintf(f,"%s",AtomValue(ttt));
		if(ListTail(li)) p+=fprintf(f,", ");
		first=0;
	}
	if(!has_GG)
	{
		if(!first) p+=fprintf(f,", ");
		fprintf(f,"GG");
	}
	fprintf(f,"\n\n");
	fclose(f);
	
	NoQuotes=1;
	
	if(OutputDirectory!=NULL)
		sprintf(cbuf,"%s/mdl_ini%d.F",OutputDirectory,fno);
	else
		sprintf(cbuf,"mdl_ini%d.F",fno);
	f=fopen(cbuf,"w");
	
	
	if(EvalPrm)
	{
		for(li=params;li;li=ListTail(li))
		{
			Term t=CompoundArg1(ListFirst(li));
			if(t && !is_integer(t) && !is_float(t))
				SetCompoundArg(ListFirst(li),1,NewFloat(EvalParameter(t)));
		}
	}
	
	fprintf(f,"*     LanHEP output produced at ");
	time(&tm);
	fprintf(f,ctime(&tm));
	if(ModelName)
		fprintf(f,"*     Model named '%s'\n",ModelName);
	
	fprintf(f,"\n      subroutine ModelConstIni(*)\n      implicit none\n\n");
	fprintf(f,"#include \"model.h\"\n\n");

	if(ExtFuncList)
	{
		fprintf(f,"      double precision ");
		for(li=ExtFuncList;li;li=ListTail(li))
			fprintf(f,"%s%c",AtomValue(ListFirst(li)),ListTail(li)?',':'\n');
		fprintf(f,"\n");
	}
	fprintf(f,"      integer m1,m2,m3,m4\n\n");				

	for(li=params;li;li=ListTail(li))
	{
		Term t, val;
		char mt[64];
		t=CompoundName(ListFirst(li));
		if(strcmp(AtomValue(t),"pi")==0)
		  continue;
		val=CompoundArg1(ListFirst(li));
		if(val==0) val=NewInteger(0);
		if(is_float(val) || is_integer(val))
		{
		  fortr_digi=1;
		  sWriteTerm(mt,val);
		  fortr_digi=0;
			fprintf(f,"      %s = %s\n",AtomValue(t),mt);
		}
	}

	
	for(li=params;li;li=ListTail(li))
	{
		Term t, val;
		char mt[1024], *mtp;
		t=CompoundName(ListFirst(li));
		if(strcmp(AtomValue(t),"pi")==0)
		  continue;
		val=CompoundArg1(ListFirst(li));
		if(val==0) val=NewInteger(0);
		if(is_float(val) || is_integer(val)) continue;
		
		p=fprintf(f,"      %s = ",AtomValue(t));
		fortr_digi=1;
		sWriteTerm(mt,val);
		fortr_digi=0;
		mtp=mt;
		while(p+strlen(mtp)>60)
		{
			char sv;
			char *sp=mtp+55-p;
			while(!isalnum(*sp)) sp++;
			while( isalnum(*sp) || *sp=='(') sp++;
			sv=(*sp);(*sp)=0;
			fprintf(f,"%s\n     &      ",mtp);
			(*sp)=sv;
			mtp=sp;
			p=12;
		}
		fprintf(f,"%s\n",mtp);
	}
	
	cls_write_matr(f);
	
	fprintf(f,"\n      end\n\n*************************************");
	fprintf(f,"**********\n\n      subroutine ModelVarIni(sqrtS, *)\n\
      implicit none\n\
      double precision sqrtS\n\
      double precision Alfas\n\
\n\
#include \"model.h\"\n\
\n\
c      double precision ALPHAS2\n\
c      external ALPHAS2\n\
\n\
c      Alfas = ALPHAS2(sqrtS)\n\
c      GG = sqrt(4*pi*Alfas)\n\
      end\n\n");
	
	fprintf(f,"************************************************\n\n");
	fprintf(f,"\
      subroutine ModelDigest\n\
      implicit none\n\
\n\
#include \"model.h\"\n\
\n");
/*	
	for(li=params;li;li=ListTail(li))
	  {
	    Term t=CompoundName(ListFirst(li));
	    fprintf(f,"      write(16,*) '%s=',%s\n",AtomValue(t),AtomValue(t));
	  } 
*/
      fprintf(f,"\n      end\n\n");
	
	fprintf(f,"      subroutine ModelDefaults\n");
    fprintf(f,"\n      end\n\n");
	
	for(li=fainclude;li;li=ListTail(li))
		fprintf(f,"#include \"%s\"\n\n",AtomValue(ListFirst(li)));
		
/*	fprintf(f,"#include \"alphas2.h\"\n\n");*/

	fclose(f);

	NoQuotes=0;
}

void FAsqparam(FILE *f)
{
	List l;
	int z=0;
	
	for(l=params;l;l=ListTail(l))
	{
		Term t=ListFirst(l);
		Term val=CompoundArg1(t);
		if(is_compound(val)&&CompoundName(val)==OPR_POW &&
				is_atom(CompoundArg1(val))&&CompoundArg2(val)==NewInteger(2))
		{
			fprintf(f,"%s^(n_?EvenQ) ^:= %s^(n/2);\n",
					AtomValue(CompoundArg1(val)),AtomValue(CompoundName(t)));
			z++;
		}
	}
	fprintf(f,"\n");	
/*	if(z)
	{
		fprintf(f,"FormSubst = \"\\\n");
		for(l=params;l;l=ListTail(l))
		{
			Term t=ListFirst(l);
			Term val=CompoundArg1(t);
			if(is_compound(val)&&CompoundName(val)==OPR_POW &&
					is_atom(CompoundArg1(val))&&CompoundArg2(val)==NewInteger(2))
			{
				z--;
				fprintf(f,"id %s^2 = %s;\\n\\\nid %s^-2 = %s^-1;\\n%c\n",
				AtomValue(CompoundArg1(val)),AtomValue(CompoundName(t)),
				AtomValue(CompoundArg1(val)),AtomValue(CompoundName(t)),
						z?'\\':'\"');
			}
		}
	}
	fprintf(f,"\n");
*/
}

int SecondVaFu=0;
	
void WriteParameters(int fno, char *name)
	{
	FILE *f1, *f2;
	List li;
	if(TexOutput)
		{
		texWriteParams(fno, name);
		return;
		}
	if(FAOutput)
		return;
	if(OutputDirectory!=NULL)
		sprintf(wpbuf,"%s/vars%d.mdl",OutputDirectory,fno);
	else
		sprintf(wpbuf,"vars%d.mdl",fno);
	f1=fopen(wpbuf,"w");
	if(f1==NULL)
		{
		printf("Can not open file \'%s\' for writing.\n",wpbuf);
		perror("");
		return;
		}
		
	if(OutputDirectory!=NULL)
		sprintf(wpbuf,"%s/func%d.mdl",OutputDirectory,fno);
	else
		sprintf(wpbuf,"func%d.mdl",fno);
	f2=fopen(wpbuf,"w");
	if(f2==NULL)
		{
		printf("Can not open file \'%s\' for writing.\n",wpbuf);
		perror("");
		return;
		}
	fprintf(f1,"%s\n Variables \n",name);
	fprintf(f2,"%s\n Constraints \n",name);
	if(MicroOmega)
		fprintf(f1," "),fprintf(f2," ");
	fprintf(f1," Name | Value       |>  Comment                                   <|\n");
	fprintf(f2," Name |> Expression");
	WriteBlank(f2,P_D_WIDTH-13);
	fprintf(f2,"<|> Comment                      <|\n");
	li=params;
	while(!is_empty_list(li))
		{
		Term ttt, val, comm;
		ttt=ListFirst(li);
		val=CompoundArg1(ttt);
		comm=CompoundArg2(ttt);
		if(val==0) val=NewInteger(0);
		if(SecondVaFu && !(is_integer(val)||is_float(val))
				&& micro_allpr(CompoundName(ttt),val))
			val=NewFloat(EvalParameter(CompoundName(ttt)));
		if(MicroOmega)
			fprintf((is_integer(val)||is_float(val))?f1:f2,
			"%c",(micro_allpr(CompoundName(ttt),val)&&SecondVaFu==0)?'*':' ');
		if(is_integer(val) || is_float(val) || EvalPrm)
			{
			int sp;
			sp=fprintf(f1,AtomValue(FunctorName(CompoundFunctor(ttt))));
			WriteBlank(f1,6-sp);
			fprintf(f1,"|");
			if(is_integer(val) || is_float(val))
				sp=fWriteTerm(f1,val);
			else
				sp=fprintf(f1,"%g",EvalParameter(CompoundName(ttt)));
			WriteBlank(f1,13-sp);
			fprintf(f1,"|");
			if(comm!=0)
				fprintf(f1,"%s\n",AtomValue(comm));
			else
				fprintf(f1,"\n");
/*			if(comm)
				sp=fprintf(f1,"%s",AtomValue(comm));
			else
				sp=0;
			WriteBlank(f1,23-sp);fprintf(f1,"|\n");*/
			}
		else
			{
			int sp;
			sp=fprintf(f2,AtomValue(FunctorName(CompoundFunctor(ttt))));
			WriteBlank(f2,6-sp);
			fprintf(f2,"|");
			NoQuotes=1;
			val=CopyTerm(val);
			if(ChepVersion>3)
				rpl_pow(val);
			sp=fWriteTerm(f2,val);
			NoQuotes=0;
			FreeAtomic(val);
			WriteBlank(f2,P_D_WIDTH-sp);
			if(sp>longest_pdline)
				longest_pdline=sp;
			fprintf(f2,"|");
			if(comm!=0)
				fprintf(f2,"%s\n",AtomValue(comm));
			else
				fprintf(f2,"\n");
/*			if(comm)
				sp=fprintf(f2,"%s",AtomValue(comm));
			else
				sp=0;
			WriteBlank(f2,47-sp);fprintf(f2,"|\n");*/
			}
		li=ListTail(li);
		}
	fclose(f1);
	fclose(f2);
	}
	
void ClearParameter(Atom p)
	{
	List l;
	l=params;
	while(!is_empty_list(l))
		{
		if(CompoundName(ListFirst(l))==p)
			{
			params=CutFromList(params,l);
			break;
			}
		l=ListTail(l);
		}
	RemoveAtomProperty(p,PROP_TYPE);
	}


	
void ChangeParameterValue(Atom t, double val)
	{
	List l;
	for(l=params;l;l=ListTail(l))
		{
		if(CompoundName(ListFirst(l))==t)
			{
			SetCompoundArg(ListFirst(l),1,NewFloat(val));
			return;
			}
		}
	puts("Internal error (prms01)");
	}	

	
double sort4(double m1, double m2, double m3, double m4, double dn)
{
	int n;
	int i,f=0;
	double m[4];
	n=(int)floor(dn-0.5);
	m[0]=m1;
	m[1]=m2;
	m[2]=m3;
	m[3]=m4;
	
	do
	{
		f=0;
		for(i=0;i<3;i++)
			if(fabs(m[i])>fabs(m[i+1]))
			{
				double tmp;
				tmp=m[i];
				m[i]=m[i+1];
				m[i+1]=tmp;
				f=1;
			}
			
	} while(f);
	
	return m[n];
}
	


static double eval_ef(Atom, int, double *);

double EvalParameter(Term t)
	{

		
	if(t==0)
		return 1.0;
		
	if(t==A_SQRT2)
		return sqrt((double)2.0);
	if(t==A_I)
		return 1.0;
	
	if(is_integer(t))
		return (double)IntegerValue(t);
	if(is_float(t))
		return FloatValue(t);
	if(is_atom(t) && is_parameter(t))
		{
		List l;
		
		
		l=GetAtomProperty(t,A_DUMMY_PRM);
		if(l)
		{
			double ret=0;
			List l1;
			for(l1=l;l1;l1=ListTail(l1))
				ret+=IntegerValue(CompoundArg1(ListFirst(l1)))*
						EvalParameter(CompoundArg2(ListFirst(l1)));
			return ret;
		}
		
		for(l=params;l;l=ListTail(l))
			{
			if(CompoundName(ListFirst(l))==t)
				{
				Term t1;
				double dd;
				t1=CompoundArg1(ListFirst(l));
				if(is_integer(t1))
					return (double)IntegerValue(t1);
				if(is_float(t1))
					return FloatValue(t1);
				if(t1==0)
					return 0.0;
				if(CompoundArgN(ListFirst(l),3)!=0)
					return FloatValue(CompoundArgN(ListFirst(l),3));
				dd=EvalParameter(t1);
				SetCompoundArg(ListFirst(l),3,NewFloat(dd));
				return dd;
				}
			}
		printf("Internal error (prm01'%s')\n",AtomValue(t));
		return 0.0;
		}	
		
	if(is_compound(t) && GetAtomProperty(CompoundName(t),A_EXT_FUNC) &&
		GetAtomProperty(CompoundName(t),CompoundName(t)))
	{
		double arr[30];
		int i;
		for(i=1;i<=CompoundArity(t);i++)
			arr[i-1]=EvalParameter(CompoundArgN(t,i));
		return eval_ef(CompoundName(t),CompoundArity(t),(double *)arr);
	}

	if(is_compound(t) && strcmp(AtomValue(CompoundName(t)),"sort4")==0
			&& CompoundArity(t)==5)
	{
		double t1,t2,t3,t4,t5;
		t1=EvalParameter(CompoundArgN(t,1));
		t2=EvalParameter(CompoundArgN(t,2));
		t3=EvalParameter(CompoundArgN(t,3));
		t4=EvalParameter(CompoundArgN(t,4));
		t5=EvalParameter(CompoundArgN(t,5));		
		t1=sort4(t1,t2,t3,t4,t5);
		return t1;
		
	}
	

	if(is_compound(t) && GetAtomProperty(CompoundName(t),A_EXT_FUNC))
	{
		return 0.0;
	}
	
	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==OPR_PLUS)
		return EvalParameter(CompoundArg1(t))+EvalParameter(CompoundArg2(t));
	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==OPR_MINUS)
		return EvalParameter(CompoundArg1(t))-EvalParameter(CompoundArg2(t));
	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==OPR_MLT)
		return EvalParameter(CompoundArg1(t))*EvalParameter(CompoundArg2(t));
	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==OPR_DIV)
		return EvalParameter(CompoundArg1(t))/EvalParameter(CompoundArg2(t));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==OPR_MINUS)
		return -EvalParameter(CompoundArg1(t));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==A_SQRT)
		return sqrt(EvalParameter(CompoundArg1(t)));
	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==A_CURT)
		return cu_rt(EvalParameter(CompoundArg1(t)),EvalParameter(CompoundArg2(t)));
	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==OPR_POW)
		return pow(EvalParameter(CompoundArg1(t)),EvalParameter(CompoundArg2(t)));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==A_SIN)
		return sin(EvalParameter(CompoundArg1(t)));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==A_COS)
		return cos(EvalParameter(CompoundArg1(t)));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==A_FABS)
		return fabs(EvalParameter(CompoundArg1(t)));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==A_TAN)
		return tan(EvalParameter(CompoundArg1(t)));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==A_ASIN)
		return asin(EvalParameter(CompoundArg1(t)));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==A_ACOS)
		return acos(EvalParameter(CompoundArg1(t)));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==A_ATAN)
		return atan(EvalParameter(CompoundArg1(t)));
	if(is_compound(t) && CompoundArity(t)==2 && CompoundName(t)==A_ATAN2)
		return atan2(EvalParameter(CompoundArg1(t)),EvalParameter(CompoundArg2(t)));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==A_LOG)
		return log(EvalParameter(CompoundArg1(t)));
	if(is_compound(t) && CompoundArity(t)==1 && CompoundName(t)==A_EXP)
		return exp(EvalParameter(CompoundArg1(t)));
	if(is_compound(t) && CompoundArity(t)==3 && CompoundName(t)==A_IF)
	{
		double t1;
		t1=EvalParameter(CompoundArg1(t));
		if(t1>0)
			t1=EvalParameter(CompoundArg2(t));
		else
			t1=EvalParameter(CompoundArgN(t,3));
		return t1;
	}	
	
	if(is_list(t))
		{
		double ret=1.0;
		List l;
		for(l=t;l;l=ListTail(l))
			{
			double tpv;
			int tpw,i;
			tpv=EvalParameter(CompoundArg1(ListFirst(l)));
			tpw=IntegerValue(CompoundArg2(ListFirst(l)));
			if(tpw>0)
				for(i=0;i<tpw;i++)
					ret*=tpv;
			else
				for(i=0;i>tpw;i--)
					ret/=tpv;
			}
		return ret;
		}
	
	printf("Error: can not evaluate ");
	WriteTerm(t);
	puts("");
	
	return 0.0;
	}

Term InterfEvalParam(Term t, Term ind)
	{
	double ret;
	if(!is_compound(t) || CompoundArity(t)!=1)
		return 0;
	ret=EvalParameter(CompoundArg1(t));
	WriteTerm(CompoundArg1(t));
	printf("=%f\n",ret);
	FreeAtomic(t);
	return 0;
	
	}
	
Term ProcTailPrm(Term t, Term ind)
{
	List l1,l2;
	
	if(!is_compound(t) || CompoundArity(t)!=1 || !is_list(CompoundArg1(t)))
	{
		ErrorInfo(112);
		puts("wrong syntax in 'tail_prm' call.");
		return 0;
	}
	
	l1=ConsumeCompoundArg(t,1);
	FreeAtomic(t);
	t=l1;
	l2=NewList();
	for(l1=t;l1;l1=ListTail(l1))
	{
		if(!is_atom(ListFirst(l1)))
		{
			ErrorInfo(112);
			printf("tail_prm: '");
			WriteTerm(ListFirst(l1));
			puts("': only parameters expected.");
			continue;
		}
		if(!is_parameter(ListFirst(l1)))
		{
			ErrorInfo(112);
			printf("tail_prm: '");
			WriteTerm(ListFirst(l1));
			puts("': it is not a parameter.");
			continue;
		}
	}
	
	for(l1=params;l1;l1=ListTail(l1))
	{
		if(ListMember(t,CompoundName(ListFirst(l1))))
		{
			Term t1;
			t1=ListFirst(l1);
			l2=AppendLast(l2,t1);
			ChangeList(l1,0);
		}
	}
	
rpt:
	for(l1=params;l1;l1=ListTail(l1))
	{
		if(ListFirst(l1)==0)
		{
			params=CutFromList(params,l1);
			goto rpt;
		}
	}
	
	params=ConcatList(params,l2);
	
	for(l1=params;l1;l1=ListTail(l1))
	{
		if(ListMember(t,CompoundName(ListFirst(l1))))
		{
			Term t1;
			double d;
			t1=ListFirst(l1);
			if(!is_integer(CompoundArg1(t1)) && !is_float(CompoundArg1(t1)))
			{
				d=EvalParameter(CompoundName(t1));
				SetCompoundArg(t1,1,NewFloat(d));
			}
		}
	}
	
	return 0;
}

Term ProcCHEPPrm(Term t, Term ind)
{
	int mono;
	char fnbuf[512];
	int has_fn=0;
	FILE *fin, *fout;
	
	if(!is_compound(t) || CompoundArity(t)<1)
	{
		ErrorInfo(113);
		puts("wrong syntax in 'read_chep_prm' statement");
		return 0;
	}
	
	if(!is_integer(CompoundArg1(t)))
	{
		ErrorInfo(113);
		puts("'read_chep_prm': the first argument should be integer.");
		return 0;
	}
	
	mono=IntegerValue(CompoundArg1(t));

	if(CompoundArity(t)>1)
		has_fn=1;
	
	if(has_fn && !is_atom(CompoundArg2(t)))
	{
		ErrorInfo(113);
		puts("'read_chep_prm': the second argument (directory) should be string const");
		return 0;
	}
	
	if(has_fn)
		sprintf(fnbuf,"%s/prm_from_chep.mdl",AtomValue(CompoundArg2(t)));
	else
		sprintf(fnbuf,"prm_from_chep.mdl");
	
	fout=fopen(fnbuf,"w");
	if(fout==NULL)
	{
		ErrorInfo(113);
		printf("'read_chep_prm': can not open output file '%s'\n",fnbuf);
		perror("'read_chep_prm':");
		return 0;
	}
	
	if(has_fn)
		sprintf(fnbuf,"%s/vars%d.mdl",AtomValue(CompoundArg2(t)),mono);
	else
		sprintf(fnbuf,"vars%d.mdl",mono);
	
	fin=fopen(fnbuf,"r");
	if(fin==NULL)
	{
		ErrorInfo(113);
		printf("'read_chep_prm': can not open input file '%s'\n",fnbuf);
		perror("'read_chep_prm':");
		return 0;
	}
	
	fprintf(fout,"%%\n%% %s\n%%\n\n",fnbuf);
	
	fgets(fnbuf,510,fin);
	fgets(fnbuf,510,fin);
	fgets(fnbuf,510,fin);
	
	while(fgets(fnbuf,510,fin)!=NULL)
	{
		int i, pos1=0, pos2=0;
		for(i=0;i<500;i++)
			if(fnbuf[i]=='|')
				break;
		
		if(i==500)
		{
			ErrorInfo(113);
			puts("'read_chep_prm': corrupted 'vars' file.");
			return 0;
		}
		
		pos1=i;
		
		for(i=pos1+1;i<500;i++)
			if(fnbuf[i]=='|')
				break;
		
		if(i==500)
		{
			ErrorInfo(113);
			puts("'read_chep_prm': corrupted 'vars' file.");
			return 0;
		}
		
		pos2=i;
		
		fnbuf[pos1]=0;
		fnbuf[pos2]=0;
		
		while(fnbuf[pos2-1]==' ')
		{
			pos2--;
			fnbuf[pos2]=0;
		}
				
		fprintf(fout,"parameter %s = %s.\n",fnbuf,fnbuf+pos1+1);
	}
	
	fclose(fin);
	
	if(has_fn)
		sprintf(fnbuf,"%s/func%d.mdl",AtomValue(CompoundArg2(t)),mono);
	else
		sprintf(fnbuf,"func%d.mdl",mono);
	
	fin=fopen(fnbuf,"r");
	if(fin==NULL)
	{
		ErrorInfo(113);
		printf("'read_chep_prm': can not open input file '%s'\n",fnbuf);
		perror("'read_chep_prm':");
		return 0;
	}
	
	fprintf(fout,"\n%%\n%% %s\n%%\n\n",fnbuf);
	
	fgets(fnbuf,510,fin);
	fgets(fnbuf,510,fin);
	fgets(fnbuf,510,fin);
	
	while(fgets(fnbuf,510,fin)!=NULL)
	{
		int i, pos1=0, pos2=0;
		for(i=0;i<500;i++)
			if(fnbuf[i]=='|')
				break;
		
		if(i==500)
		{
			ErrorInfo(113);
			puts("'read_chep_prm': corrupted 'vars' file.");
			return 0;
		}
		
		pos1=i;
		
		for(i=pos1+1;i<500;i++)
			if(fnbuf[i]=='|')
				break;
		
		if(i==500)
		{
			ErrorInfo(113);
			puts("'read_chep_prm': corrupted 'vars' file.");
			return 0;
		}
		
		pos2=i;
		
		fnbuf[pos1]=0;
		fnbuf[pos2]=0;
		
		while(fnbuf[pos2-1]==' ')
		{
			pos2--;
			fnbuf[pos2]=0;
		}
		
		fprintf(fout,"parameter %s = %s.\n",fnbuf,fnbuf+pos1+1);
	}
	
	fclose(fin);
	fprintf(fout,"\n");
	fclose(fout);
	
	if(has_fn)
		sprintf(fnbuf,"%s/prm_from_chep.mdl",AtomValue(CompoundArg2(t)));
	else
		sprintf(fnbuf,"prm_from_chep.mdl");
	
	ReadFile(fnbuf);
	
	FreeAtomic(t);
	
	return 0;
}

static struct efel
	{
	int no;
	Atom fname, ffile;
	void *handler;
	extfunc func;
	struct efel *next;
	} *flist=NULL;

static int extfuncno=0;

Term ProcExtFunc(Term t, Term ind)
{
	Atom fl=0;
	if(!is_compound(t)|| CompoundArity(t)<2 || CompoundArity(t)>3 || 
		!is_atom(CompoundArg1(t))|| !is_integer(CompoundArg2(t)) ||
		(CompoundArity(t)==3 && !is_atom(CompoundArgN(t,3))))
	{
		ErrorInfo(228);
		printf("Illegal syntax in external_func statement.\n");
		return 0;
	}
	if(IntegerValue(CompoundArg2(t))<0 || IntegerValue(CompoundArg2(t))>29)
	{
		ErrorInfo(229);
		puts("external_func: numer of arguments must be from 1 to 29.");
		return 0;
	}

	if(CompoundArity(t)==3)
	{
		struct efel *ep;
		void *ha;
		extfunc hf;
		fl=CompoundArgN(t,3);
		for(ep=flist;ep;ep=ep->next)
			if(ep->ffile==fl)
				break;
		if(ep==0)
		{
			ha=dlopen(AtomValue(fl),RTLD_LAZY);
			if(ha==0)
				{
				char cbuf[128];
				sprintf(cbuf,"./%s",AtomValue(fl));
				ha=dlopen(cbuf,RTLD_LAZY);
				/*if(ha==0){dlerror();ha=dlopen(AtomValue(fl),RTLD_LAZY);}*/
				}
			if(ha==0)
				{
				ErrorInfo(9);
				printf("external_func: %s\n",dlerror());
				fl=0;
				goto cnt;
				}
		}
		else
			ha=ep->handler;
		dlerror();
		hf=(extfunc)dlsym(ha,AtomValue(CompoundArg1(t)));
		if(hf==NULL)
		{
			ErrorInfo(10);
			printf("external_func: %s\n",dlerror());
			if(ep==0)
				dlclose(ha);
			fl=0;
			goto cnt;
		}
		extfuncno++;
		ep=(struct efel *)malloc(sizeof(struct efel));
		ep->no=extfuncno;
		ep->fname=CompoundArg1(t);
		ep->ffile=fl;
		ep->handler=ha;
		ep->func=hf;
		ep->next=flist;
		flist=ep;
	}
	
cnt:
	SetAtomProperty(CompoundArg1(t), A_EXT_FUNC, CompoundArg2(t));
	if(fl)
		SetAtomProperty(CompoundArg1(t),CompoundArg1(t),NewInteger(extfuncno));
	ExtFuncList=AppendLast(ExtFuncList,CompoundArg1(t));
	return 0;
}

static double eval_ef(Atom fname, int ar, double *a)
	{
	struct efel *e;
	extfunc f;
	for(e=flist;e;e=e->next)
		if(e->fname==fname)
			break;
/*	printf("calling %s...\n",AtomValue(fname));*/
	if(e==0)
		{
		puts("Internal error (evlef01)");
		return 0.0;
		}
	f=e->func;
	switch(ar)
		{
		case 1: return f(a[0]);
		case 2: return f(a[0],a[1]);
		case 3: return f(a[0],a[1],a[2]);
		case 4: return f(a[0],a[1],a[2],a[3]);
		case 5: return f(a[0],a[1],a[2],a[3],a[4]);
		case 6: return f(a[0],a[1],a[2],a[3],a[4],a[5]);
		case 7: return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6]);
		case 8: return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7]);
		case 9: return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8]);
		case 10:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9]);
		
		case 11:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10]);
		case 12:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11]);
		case 13:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12]);
		case 14:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13]);
		case 15:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14]);
		case 16:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15]);
		case 17:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16]);
		case 18:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17]);
		case 19:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18]);		
		case 20:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19]);
		case 21:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
				a[20]);
		case 22:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
				a[20],a[21]);
		case 23:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
				a[20],a[21],a[22]);
		case 24:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
				a[20],a[21],a[22],a[23]);
		case 25:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
				a[20],a[21],a[22],a[23],a[24]);
		case 26:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
				a[20],a[21],a[22],a[23],a[24],a[25]);
		case 27:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
				a[20],a[21],a[22],a[23],a[24],a[25],a[26]);
		case 28:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
				a[20],a[21],a[22],a[23],a[24],a[25],a[26],a[27]);
		case 29:return f(a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],
				a[10],a[11],a[12],a[13],a[14],a[15],a[16],a[17],a[18],a[19],
				a[20],a[21],a[22],a[23],a[24],a[25],a[26],a[27],a[28]);
		default:
				puts("Internal error (evlef02)");
				return 0.0;
		}
	}	




Term rm_zero(Term t)
{
	int i;
	if(!is_compound(t))
		return t;
	
	return t;
	
	for(i=1;i<=CompoundArity(t);i++)
	{
		SetCompoundArg(t,i,rm_zero(ConsumeCompoundArg(t,i)));
	}
	
	if(CompoundName(t)==OPR_POW && CompoundArity(t)==2 && 
			CompoundArg1(t)==NewInteger(0))
	{
		FreeAtomic(t); 
		return NewInteger(0);
	}
	
	if(CompoundName(t)==OPR_MLT && CompoundArity(t)==2 && 
			(CompoundArg1(t)==NewInteger(0) || CompoundArg2(t)==NewInteger(0)))
	{
		FreeAtomic(t);
		return NewInteger(0);
	}
	if(CompoundName(t)==OPR_PLUS && CompoundArity(t)==2 &&CompoundArg2(t)==NewInteger(0))
	{
		Term t1=ConsumeCompoundArg(t,1);
		FreeAtomic(t);
		return rm_zero(t1);
	}
	if(CompoundName(t)==OPR_PLUS && CompoundArity(t)==2 &&CompoundArg1(t)==NewInteger(0))
	{
		Term t1=ConsumeCompoundArg(t,2);
		FreeAtomic(t);
		return rm_zero(t1);
	}
	if(CompoundName(t)==OPR_MINUS && CompoundArity(t)==2 &&CompoundArg2(t)==NewInteger(0))
	{
		Term t1=ConsumeCompoundArg(t,1);
		FreeAtomic(t);
		return rm_zero(t1);
	}
	if(CompoundName(t)==OPR_MINUS && CompoundArity(t)==2 &&CompoundArg1(t)==NewInteger(0))
	{
		Term t1=ConsumeCompoundArg(t,1);
		FreeAtomic(t);
		t1=rm_zero(t1);
		return t1==NewInteger(0)?NewInteger(0):MakeCompound1(OPR_MINUS,t1);
	}
	
	return t;
}

