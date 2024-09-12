#include <stdio.h>
#include <string.h>
#include "lanhep.h"

static Term dif_term(Term t)
	{
	if(t==A_SQRT2 || t==A_I || is_float(t) || is_integer(t))
		return 0;
	
	if(is_atom(t) && is_parameter(t))
		return MakeCompound1(OPR_USCORE,t);
	
	if(!is_compound(t))
		{
		printf("Can't variate "); WriteTerm(t); puts(""); return 0;
		}
	
	if(CompoundName(t)==OPR_PLUS && CompoundArity(t)==2)
		{
		Term d1,d2;
		d1=dif_term(CompoundArg1(t));
		d2=dif_term(CompoundArg2(t));
		if(d1==0 && d2==0)
			return 0;
		if(d1==0)
			return d2;
		if(d2==0)
			return d1;
		return MakeCompound2(OPR_PLUS,d1,d2);
		}
	if(CompoundName(t)==OPR_MINUS && CompoundArity(t)==2)
		{
		Term d1,d2;
		d1=dif_term(CompoundArg1(t));
		d2=dif_term(CompoundArg2(t));
		if(d1==0 && d2==0)
			return 0;
		if(d1==0)
			return d2;
		if(d2==0)
			return d1;
		return MakeCompound2(OPR_MINUS,d1,d2);
		}
	if(CompoundName(t)==OPR_MINUS && CompoundArity(t)==1)
		{
		Term d1;
		d1=dif_term(CompoundArg1(t));
		if(d1==0)
			return 0;
		return MakeCompound1(OPR_MINUS,d1);
		}
		
	if(CompoundName(t)==A_SQRT && CompoundArity(t)==1)
		{
		Term d1;
		d1=dif_term(CompoundArg1(t));
		if(d1==0)
			return 0;
		if(is_compound(d1) && CompoundName(d1)==OPR_MLT &&
			CompoundArg1(d1)==NewInteger(2))
				{
				Term d2;
				d2=ConsumeCompoundArg(d1,2);
				FreeAtomic(d1);
				return MakeCompound1(OPR_MINUS, 
					MakeCompound2(OPR_DIV,d2,t));
				}
		return MakeCompound1(OPR_MINUS,
			MakeCompound2(OPR_DIV, d1, 
				MakeCompound2(OPR_MLT, NewInteger(2), t)));
		}
		
	if(CompoundName(t)==OPR_POW && CompoundArity(t)==2)
		{
		Term d1;
		int ppow;
		ppow=IntegerValue(CompoundArg2(t));
		d1=dif_term(CompoundArg1(t));
		if(d1==0)
			return 0;
		if(ppow==2)
			return MakeCompound2(OPR_MLT,NewInteger(2),
				MakeCompound2(OPR_MLT, CompoundArg1(t), d1));
		return MakeCompound2(OPR_MLT, d1, MakeCompound2(OPR_POW,
			CompoundArg1(t),NewInteger(ppow-1)));
		}
	
	if(CompoundName(t)==OPR_MLT && CompoundArity(t)==2)
		{
		Term d1,d2;
		d1=dif_term(CompoundArg1(t));
		d2=dif_term(CompoundArg2(t));
		if(d1==0 && d2==0)
			return 0;
		if(d1==0)
			return MakeCompound2(OPR_MLT,CompoundArg1(t),d2);
		if(d2==0)
			return MakeCompound2(OPR_MLT,CompoundArg2(t),d1);
		return MakeCompound2(OPR_PLUS,
			MakeCompound2(OPR_MLT,CompoundArg1(t),d2),
			MakeCompound2(OPR_MLT,CompoundArg2(t),d1));
		}
		
	if(CompoundName(t)==OPR_DIV && CompoundArity(t)==2)
		{
		Term d1,d2;
		d1=dif_term(CompoundArg1(t));
		d2=dif_term(CompoundArg2(t));
		if(d1==0 && d2==0)
			return 0;
		if(d1==0)
			return 
				MakeCompound1(OPR_MINUS,
					MakeCompound2(OPR_DIV,
						MakeCompound2(OPR_MLT,d2,CompoundArg1(t)),
						MakeCompound2(OPR_POW,CompoundArg2(t),NewInteger(2))
							     )
							 );
		if(d2==0)
			return MakeCompound2(OPR_DIV,d1,CompoundArg2(t));
				
		return 
			MakeCompound2(OPR_DIV,
				MakeCompound2(OPR_MINUS,
					MakeCompound2(OPR_MLT,d1,CompoundArg2(t)),
					MakeCompound2(OPR_MLT,d2,CompoundArg1(t))
							 ),
				MakeCompound2(OPR_POW,CompoundArg2(t),NewInteger(2))
				);
		}
	
	printf("Can't variate "); WriteTerm(t); puts("");
	return 0;
	
	}

Term VarVer(Term t, Term ind)
	{
	List l,l1,l2;
	puts("Parameters");
	l=all_param_list();
	for(l1=l;!is_empty_list(l1);l1=ListTail(l1))
		{
		Term aa;
		aa=CompoundArg1(ListFirst(l1));
		if(!is_float(aa) && !is_integer(aa))
			{
			printf("%s: ",AtomValue(CompoundName(ListFirst(l1))));
			WriteTerm(aa);
			printf(" ->  ");
			WriteTerm(dif_term(aa));
			puts("");
			}
		}
		
	puts("\nVertices");

	l2=l=all_vert_list();
		
	for(l1=l;!is_empty_list(l1);l1=ListTail(l1))
		{
		Term a2;
		List l,lp,lm;
		a2=CopyTerm(ListFirst(l1));
		alg2_symmetrize(a2);
		alg2_common_s(a2);
		alg2_common_n(a2);
		alg2_red_cos(a2);
		alg2_red_orth(a2);
		if(CompoundArg2(a2)==NewInteger(0) ||
			 is_empty_list(CompoundArgN(a2,5)))
				continue;


		WriteVertex(CompoundArg1(a2));
		printf(" ");
		if(CompoundArgN(a2,3))
			WriteTerm(CompoundArgN(a2,3));
		else
			printf("[]");
		printf(" ");
		lp=lm=NewList();
		for(l=CompoundArgN(a2,3);!is_empty_list(l);l=ListTail(l))
			{
			Term aa;
			aa=ListFirst(l);
			if(IntegerValue(CompoundArg2(aa))<0)
				{
				aa=CopyTerm(aa);
				SetCompoundArg(aa,2,NewInteger(-IntegerValue(CompoundArg2(aa))));
				if(IntegerValue(CompoundArg2(aa))==1)
					lm=AppendLast(lm,CompoundArg1(aa));
				else
					lm=AppendLast(lm,aa);
				}
			else
				{
				if(IntegerValue(CompoundArg2(aa))==1)
					lp=AppendLast(lp,CompoundArg1(aa));
				else
					lp=AppendLast(lp,aa);
				}
			}
			
		if(lp)
			lp=l2mult(lp);
		if(lm)
			lm=l2mult(lm);
		if(lp==0 && lm==0)
			{
			printf("0");
			goto cnt;
			}
		if(lm==0)
			{
			WriteTerm(dif_term(lp));
			goto cnt;
			}
		if(lp==0)
			{
			WriteTerm(dif_term(MakeCompound2(OPR_DIV,NewInteger(1),lm)));
			goto cnt;
			}
		WriteTerm(dif_term(MakeCompound2(OPR_DIV,lp,lm)));
		
		cnt:
		puts("");
		FreeAtomic(a2);
		}
	FreeAtomic(l2);
	return 0;
	}

extern int LagrHashSize;
extern List *lagr_hash;

Term ProcChVertex(Term t, List ind)
{
	List l, pl, ml;
	Term rpl;
	Atom a1, a2;
	int pli;
	int ii;
	
	if(lagr_hash==NULL)
	{
		ErrorInfo(107);
		puts("ChVertex: no vertices");
		return 0;
	}
	
	
	if(!is_compound(t)||CompoundArity(t)!=2)
	{
		ErrorInfo(107);
		puts("wrong call to ChVertex");
		return 0;
	}
	
	pl=CompoundArg1(t);
	rpl=CompoundArg2(t);
	if(!is_list(pl)|| !is_compound(rpl) || CompoundArity(rpl)!=2)
	{
		ErrorInfo(107);
		puts("wrong call to ChVertex");
		return 0;
	}
	a1=CompoundArg1(rpl);a2=CompoundArg2(rpl);
	if(!is_parameter(a1)||!is_parameter(a2))
	{
		ErrorInfo(107);
		puts("wrong call to ChVertex");
		return 0;
	}
	
	pli=ListLength(pl);
	
	for(ii=0;ii<LagrHashSize;ii++)
	{
	
		for(l=lagr_hash[ii];l;l=ListTail(l))
		{
			List cpl=CompoundArg1(ListFirst(l));
			List ll1,ll2;

			if(cpl==0 || ListLength(cpl)!=pli)
				continue;

			for(ll1=pl,ll2=cpl;ll1;ll1=ListTail(ll1),ll2=ListTail(ll2))
				if(ListFirst(ll1)!=CompoundArg1(ListFirst(ll2)))
					break;
			if(is_empty_list(ll1))
				break;
		}
	if(l)
		break;
	}
	
	if(is_empty_list(l))
	{
	WarningInfo(108);puts("ChVertex: vertex not found");
	return 0;
	}
	
	ml=CompoundArgN(ListFirst(l),5);
	ii=0;
		
	for(l=ml;l;l=ListTail(l))
	{
		List l1;
		for(l1=CompoundArg2(ListFirst(l));l1;l1=ListTail(l1))
			if(CompoundArg1(ListFirst(l1))==a1)
			{
				SetCompoundArg(ListFirst(l1),1,a2);ii++;
			}
	}
	
	
	if(ii==0)
	{
	WarningInfo(107);printf("ChVertex: vertex has no '%s' within\n",
		AtomValue(a1));
	}	
	
	return 0;
}

static int inifile=0;

Term ProcMkProc(Term t, Term ind)
{
	Atom prt[4];
	Atom mass[4];
	Term color[4];
	int spin[4];
	int i;
	int neufact=1;
	double thcut=0.0;
	char pname[128];
	FILE *f;
	
	if(CompoundArity(t)<4)
		{
		ErrorInfo(2000);
		puts("mkProc: wrong argument number.");
		return 0;
		}
		
	for(i=5;i<=CompoundArity(t);i++)
		{
		Term t1=CompoundArgN(t,i);
		if(is_compound(t1) && is_atom(CompoundArg1(t1)) &&
			strcmp(AtomValue(CompoundArg1(t1)),"THETACUT")==0)
			{
			if(is_integer(CompoundArg2(t1)))
				thcut=IntegerValue(CompoundArg2(t1));
			else if(is_float(CompoundArg2(t1)))
				thcut=FloatValue(CompoundArg2(t1));
			else
				{
				ErrorInfo(303);puts("wrong THETACUT value.");
				continue;
				}
			continue;
			}
		ErrorInfo(304);
		puts("wrong argument in mkProc.");
		}
			
		
		
	for(i=0;i<4;i++)
		{
		Term prp, t7;
		prt[i]=CompoundArgN(t,i+1);
		if(!is_particle(prt[i],NULL))
			{
			ErrorInfo(2001);
			WriteTerm(prt[i]);
			puts(": is not a particle.");
			return 0;
			}
		prp=GetAtomProperty(prt[i],PROP_TYPE);
		if(CompoundName(prp)!=OPR_PARTICLE)
			{
			ErrorInfo(2001);
			WriteTerm(prt[i]);
			puts(": is not a particle.");
			return 0;
			}
		if(prt[i]==CompoundArg2(prp))
			prp=GetAtomProperty(CompoundArg1(prp),PROP_TYPE);
		spin[i]=IntegerValue(CompoundArgN(prp,4));
		mass[i]=CompoundArgN(prp,5);
		color[i]=GetAtomProperty(prt[i],A_COLOR);
		t7=CompoundArgN(prp,7);
		if(i<2 && (t7==A_LEFT||t7==A_RIGHT))
			neufact*=2;
		}
		
	sprintf(pname,"%s%s%s%s",AtomValue(prt[0]),AtomValue(prt[1]),
			AtomValue(prt[2]),AtomValue(prt[3]));
	
	for(i=0;pname[i];i++)
		{
		if(pname[i]=='~') pname[i]='_';
		if(pname[i]=='+') pname[i]='p';
		if(pname[i]=='-') pname[i]='m';
		}
	
	f=fopen("scan.bat",inifile?"a":"w");
	if(f==NULL)
		{
		ErrorInfo(2000);
		puts("mkProc: can not open scan.bat");
		return 0;
		}
		
	if(!inifile)
		{
		fprintf(f,"#!/bin/sh\n\n");
		inifile=1;
		}
		
	fprintf(f,"echo Generating process %s   `date`\n\n",pname);
	fprintf(f,"echo Process %s:  >> scan.log\n",pname);
	
	fprintf(f,"cat > proc.m <<END\n");
	fprintf(f,"process = {prt[\"%s\"],prt[\"%s\"]} ->",
			AtomValue(prt[0]),AtomValue(prt[1]));
	fprintf(f," {prt[\"%s\"],prt[\"%s\"]}\n",
			AtomValue(prt[2]),AtomValue(prt[3]));
	fprintf(f,"dir = SetupCodeDir[\"scan_%s\"]\n",pname);
	fprintf(f,"SetOptions[InsertFields,Model->model%d, GenericModel->model%d,\n",
				ModelNumber,ModelNumber);
	fprintf(f,"       ExcludeParticles->{ ");
	if(color[0]&&color[1]&&color[2]&&color[3])
		{
		int glu=1,gno=1;
		if(GetAtomProperty(prt[0],A_ANTI)!=prt[1]) glu=0;
		if(spin[0]==0&&spin[1]==0&&spin[2]==0&&spin[3]==0) gno=0;
		if(glu)
			fprintf(f,"prt[\"G\"]%c ",gno?',':' ');
		if(gno)
			fprintf(f,"prt[\"~%c\"] ",ModelNumber>30?'G':'g');
		}
	fprintf(f,"} ]\n");
	
	fprintf(f,"END\n\n");
	
	fprintf(f,"if test ! -d scan_%s/squaredme ;\n",pname);
	fprintf(f,"then  math < scan.m;\n");
	fprintf(f,"fi\n\n");
	
	fprintf(f,"if test ! -d scan_%s/squaredme ;\n",pname);
	fprintf(f,"then echo Output directory is not created | tee -a scan.log && exit;\n");
	fprintf(f,"fi\n\n");
	
	fprintf(f,"cat > scan_%s/process.h <<END\n",pname);
	for(i=1;i<=4;i++)
		{
		fprintf(f,"#define TYPE%d %s\n",i,
				spin[i-1]==0?"SCALAR":(spin[i-1]==1?"FERMION":
					(mass[i-1]==0?"PHOTON":"VECTOR")));
		fprintf(f,"#define MASS%d %s\n",i,mass[i-1]?AtomValue(mass[i-1]):"0");
		fprintf(f,"#define CHARGE%d 0\n\n",i);
		}
	fprintf(f,"#define IDENTICALFACTOR %s\n",(prt[2]==prt[3])?"0.5":"1");
	fprintf(f,"#define COLOURFACTOR %dD0",neufact);
	if(color[0] && color[1]) fprintf(f,"/9D0");
	else if(color[0] || color[1]) fprintf(f,"/3D0");
	fprintf(f,"\n");
	fprintf(f,"#define NCOMP 2\n#include \"2to2.F\"\nEND\n\n");
	
	fprintf(f,"cp model%d.h scan_%s/model.h\n",ModelNumber,pname);
	fprintf(f,"cp mdl_ini%d.F scan_%s/mdl_ini.F\n\n",ModelNumber,pname);
	fprintf(f,"cp main.F scan_%s/\n\n",pname);
	if(thcut!=0.0)
		{
		fprintf(f,"echo \"#define THETACUT (%f*degree)\" > scan_%s/run.F\n",
					thcut,pname);
		fprintf(f,"grep -v THETACUT drivers/run.F >> scan_%s/run.F\n",pname);
		} 
	
	fprintf(f,"cd scan_%s\n",pname);
	
	fprintf(f,"if test ! -f run ;\n");
	fprintf(f,"then ./configure && cp ../makefile . ;\n");
	fprintf(f,"fi\n\n");
	fprintf(f,"rm run ru*.01000*/*\n");
	fprintf(f,"gmake\n");
	fprintf(f,"if test ! -f run  ;\n");
	fprintf(f,"then echo Run file is not created | tee -a ../scan.log && exit;\n");
	fprintf(f,"fi\n");
	fprintf(f,"./run uuuu 1000,1000\ncd ..\n");
	fprintf(f,"grep \"|    1000.000\" scan_%s/ru*.01000*/* >> scan.log\n\n\n",
				pname);
	
	fclose(f);
	return 0;
	
	}
	
	
