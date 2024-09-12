#include <setjmp.h>
#include "lanhep.h"

extern jmp_buf alg1_jmp_buf;


static void repl_ind(Term a2, List oi, List ni)
	{

	a2=CompoundArg1(a2);
	while(!is_empty_list(a2))
		{
		List t;
		t=CompoundArgN(ListFirst(a2),3);
		while(!is_empty_list(t))
			{
			List t1;
			t1=CompoundArg1(ListFirst(t));

			while(!is_empty_list(t1))
				{
				Term a;
				List loi,lni;
				a=CompoundArg2(ListFirst(t1));
				loi=oi;
				lni=ni;
				while(!is_empty_list(loi))
					{
					if(a==CompoundArg2(ListFirst(loi)))
						{
						SetCompoundArg(ListFirst(t1),2,
								CompoundArg2(ListFirst(lni)));
						break;
						}
					loi=ListTail(loi);
					lni=ListTail(lni);
					}
				t1=ListTail(t1);
				}

			t=ListTail(t);
			}
		a2=ListTail(a2);
		}				
	return ;
	}			
		

static List mk_let(Term m1, List cut, Term a1)
	{
	List l;
	List lb,le;
	List sd;
	int num1,den1;
	
	l=ConsumeCompoundArg(a1,1);
	FreeAtomic(a1);
	a1=l;
	
	
	num1=IntegerValue(ConsumeCompoundArg(m1,1));
	den1=IntegerValue(ConsumeCompoundArg(m1,2));
	l=ConsumeCompoundArg(m1,3);
	sd=ConsumeCompoundArg(m1,4);
	lb=ListSplit(l,cut,&le);
	FreeAtomic(cut);
	FreeAtomic(m1);
	
	l=a1;
	while(!is_empty_list(l))
		{
		int n1,n2,d1,d2,num,den,cf;
		List lb1,le1,lm;
		m1=ListFirst(l);
		lm=ConsumeCompoundArg(m1,3);
		lb1=CopyTerm(lb);
		le1=CopyTerm(le);
		lb1=ConcatList(lb1,lm);
		lb1=ConcatList(lb1,le1);
		SetCompoundArg(m1,3,lb1);
		
		lm=ConsumeCompoundArg(m1,4);
		lm=ConcatList(lm,CopyTerm(sd));
		SetCompoundArg(m1,4,lm);
		
		n1=num1;
		d1=den1;
		n2=IntegerValue(CompoundArg1(m1));
		d2=IntegerValue(CompoundArg2(m1));
		num=n1*n2;
		den=d1*d2;
		if(den<0) den=-den;
		cf=gcf(num,den);
		num/=cf;
		den/=cf;
		if(d1<0 && d2<0)
			{
			num=-num;
			}
		else
			if((d1<0 && d2>0) || (d1>0 && d2<0))
				{
				den=-den;
				}
		SetCompoundArg(m1,1,NewInteger(num));
		SetCompoundArg(m1,2,NewInteger(den));
		l=ListTail(l);
		}
	FreeAtomic(lb);
	FreeAtomic(le);
	FreeAtomic(sd);
	return a1;
	}



static List mk_let_d(Term m1,  Term a1)
	{
	List l;
	List lb;
	List sd;
	int num1,den1;


	l=ConsumeCompoundArg(a1,1);
	FreeAtomic(a1);
	a1=l;
	
	num1=IntegerValue(ConsumeCompoundArg(m1,1));
	den1=IntegerValue(ConsumeCompoundArg(m1,2));
	l=ConsumeCompoundArg(m1,3);
	sd=ConsumeCompoundArg(m1,4);
	lb=l;
	FreeAtomic(m1);
	
	l=a1;
	while(!is_empty_list(l))
		{
		int n1,n2,d1,d2,num,den,cf;
		List lb1,lm;
		m1=ListFirst(l);
		lm=ConsumeCompoundArg(m1,3);
		lb1=CopyTerm(lb);
		lb1=ConcatList(lb1,lm);
		SetCompoundArg(m1,3,lb1);
		
		lm=ConsumeCompoundArg(m1,4);
		lm=ConcatList(lm,CopyTerm(sd));
		SetCompoundArg(m1,4,lm);
		
		n1=num1;
		d1=den1;
		n2=IntegerValue(CompoundArg1(m1));
		d2=IntegerValue(CompoundArg2(m1));
		num=n1*n2;
		den=d1*d2;
		if(den<0) den=-den;
		cf=gcf(num,den);
		num/=cf;
		den/=cf;
		if(d1<0 && d2<0)
			{
			num=-num;
			}
		else
			if((d1<0 && d2>0) || (d1>0 && d2<0))
				{
				den=-den;
				}
		SetCompoundArg(m1,1,NewInteger(num));
		SetCompoundArg(m1,2,NewInteger(den));
		l=ListTail(l);
		}
	FreeAtomic(lb);
	FreeAtomic(sd);
	return a1;
	}


static List s_l_1(Term m1)
	{
	List l,l1;
	
	l=CompoundArgN(m1,3);

	while(!is_empty_list(l))
		{
		Term t1;
		t1=ListFirst(l);
		if(CompoundName(t1)==A_ALG1)
			{
			Term sub,a1,ila,ill;
			
			a1=sub=ConsumeCompoundArg(t1,2);
			ill=ConsumeCompoundArg(t1,1);
			FreeAtomic(t1);
			ChangeList(l,0);
					
			ila=CopyTerm(CompoundArg2(sub));
			
			repl_ind(a1,ila,ill);

			return SetIntAlgs(mk_let(m1,l,a1));
			
			}
		l=ListTail(l);
		}
	
	
	l=CompoundArgN(m1,4);
	while(!is_empty_list(l))
		{
		Term t1;
		t1=ListFirst(l);
		if(CompoundName(t1)==A_ALG1)
			{
			Term sub, a1;
			sub=alg1_inv_alg(CompoundArg2(t1));
			if(sub==0)
			{
				ErrorInfo(375);
				puts("polynomial in denominator");
				longjmp(alg1_jmp_buf,1);
			}
			
			l1=ConsumeCompoundArg(m1,4);
			l1=CutFromList(l1,l);
			SetCompoundArg(m1,4,l1);
			
			a1=sub;	
			
			return SetIntAlgs(mk_let_d(m1,a1));
			
			}
		l=ListTail(l);
		}
		
	

	
	return AppendFirst(NewList(),m1);
	}

static Term alg1_df_df(Term m, int pos, List il)
{
	List nl,l1,l2;
	Term obj;
	
	nl=ConsumeCompoundArg(m,3);
	
	l1=ListNthList(nl,pos);
	
	obj=ListFirst(l1);
	ChangeList(l1,0);
	nl=CutFromList(nl,l1);
	
	l1=ConsumeCompoundArg(obj,1);
	FreeAtomic(obj);
	obj=l1;
	
	for(l1=obj,l2=il;l1;l1=ListTail(l1),l2=ListTail(l2))
		nl=AppendLast(nl,MakeCompound2(OPR_SPECIAL,
				MakeList2(ListFirst(l1),CopyTerm(ListFirst(l2))),A_DELTA));
	
	SetCompoundArg(m,3,nl);
	RemoveList(obj);
	
	return m;
}
	
			
static Term alg1_df(Term a1, Atom p)
{
	Term ind1, ind2,l1,l2,mo,mn=0;
	ind1=ConsumeCompoundArg(a1,2);
	ind2=CopyTerm(GetAtomProperty(p,PROP_INDEX));
	
	for(l1=ind2;l1;l1=ListTail(l1))
	{
		Label l;
		Term g, r1,r2;
		l=NewLabel();
		SetCompoundArg(ListFirst(l1),2,l);
		g=CompoundArg1(ListFirst(l1));
		r1=ConsumeCompoundArg(g,1);
		r2=ConsumeCompoundArg(g,2);
		SetCompoundArg(g,1,r2);
		SetCompoundArg(g,2,r1);
	}
	
	mo=ConsumeCompoundArg(a1,1);
	
	for(l1=mo;l1;l1=ListTail(l1))
	{
		Term mt;
		int i;
		mt=ListFirst(l1);
		for(l2=CompoundArgN(mt,4);l2;l2=ListTail(l2))
			if(CompoundArg2(ListFirst(l2))==p)
			{
				ErrorInfo(371);
				printf("df: object '%s' in denominator.\n",AtomValue(p));
				FreeAtomic(mo);
				return a1;
			}
		
		for(l2=CompoundArgN(mt,3),i=1;l2;l2=ListTail(l2),i++)
			if(CompoundArg2(ListFirst(l2))==p)
				mn=AppendLast(mn,alg1_df_df(CopyTerm(mt),i,ind2));
				
		
	}
	
	FreeAtomic(mo);
	
	SetCompoundArg(a1,1,mn);
	SetCompoundArg(a1,2,ConcatList(ind1,ind2));

	return a1;
}			

void alg1_kl_to_ia(List al)
{
	for(;al;al=ListTail(al))
	{
		List l1;
		for(l1=CompoundArgN(ListFirst(al),3);l1;l1=ListTail(l1))
		{
			Term prp;
			if(CompoundName(ListFirst(l1))==OPR_LET &&
				(prp=GetAtomProperty(CompoundArg2(ListFirst(l1)),A_KEEP_LETS))
				&& CompoundArg1(prp))
			{
				List il,kll;
				kll=CopyTerm(CompoundArg1(prp));
				il=ConsumeCompoundArg(ListFirst(l1),1);
				FreeAtomic(ListFirst(l1));
				ChangeList(l1,MakeCompound2(A_ALG1,il,kll));
			}
		}
	}
}


List SetIntAlgs(List l)
	{
	List l1,lr;
	
	
	lr=NewList();
	l1=l;
	while(!is_empty_list(l1))
		{
		lr=ConcatList(lr,s_l_1(ListFirst(l1)));
		l1=ListTail(l1);
		}
	RemoveList(l);	
	return lr;
	}
	
	

Term ProcDF(Term t, Term ind)
{
	
	Term t1,t2,t3;
	int ano;
	
	if(!is_compound(t) || CompoundArity(t)<2 || CompoundArity(t)>3 ||
			!is_atom(CompoundArg1(t)) || !is_atom(CompoundArg2(t)) ||
			(CompoundArity(t)==3 && !is_atom(CompoundArgN(t,3))))
	{
		ErrorInfo(367);
		printf("wrong arguments in 'df' call.\n");
		return 0;
	}
	
	ano=CompoundArity(t)-1;
	
	t1=CompoundArg1(t);
	t2=CompoundArg2(t);
	if(ano==2)
		t3=CompoundArgN(t,3);
	else
		t3=0;
	
	FreeAtomic(t);
	
	if(GetAtomProperty(t1,PROP_TYPE)==0)
	{
		ErrorInfo(368);
		printf("undefined object '%s'\n",AtomValue(t1));
		return 0;
	}
	
	if(GetAtomProperty(t2,PROP_TYPE)==0)
	{
		ErrorInfo(368);
		printf("undefined object '%s'\n",AtomValue(t2));
		return 0;
	}
	
	if(t3 && GetAtomProperty(t3,PROP_TYPE)==0)
	{
		ErrorInfo(368);
		printf("undefined object '%s'\n",AtomValue(t3));
		return 0;
	}
	
	t=GetAtomProperty(t1,A_KEEP_LETS);
	if(t==0 || CompoundArg1(t)==0)
	{
		ErrorInfo(369);
		printf("'df': object '%s' is not defined.\n",AtomValue(t1));
		return 0;
	}
	
	if(!GetAtomProperty(t2,PROP_TYPE))
	{
		ErrorInfo(369);
		printf("'df': object '%s' is not defined.\n",AtomValue(t2));
		return 0;
	}
	
	if(t3 && !GetAtomProperty(t3,PROP_TYPE))
	{
		ErrorInfo(369);
		printf("'df': object '%s' is not defined.\n",AtomValue(t3));
		return 0;
	}
	
	t=CopyTerm(CompoundArg1(t));
	
	t=alg1_df(t,t2);
	if(t3)
		t=alg1_df(t,t3);
	
	
	if(CompoundArg1(t)==0)
	{
		FreeAtomic(t);
		return NewInteger(0);
	}
	else
	{
		if(ind)
			t=il_to_caret(t,ind);
		
		return t;
	}
	
}


