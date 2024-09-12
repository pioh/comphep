#include <setjmp.h>
#include "lanhep.h"

extern jmp_buf alg1_jmp_buf;


Term namely_w_t(Term t)
	{
	Term tt,sub;
	List l1,l2,l3;
	l3=NewList();
	tt=ConsumeCompoundArg(t,1);
	sub=ConsumeCompoundArg(t,2);
	FreeAtomic(t);
	l1=OperToList(sub,OPR_SECO);
	l2=l1;
	while(!is_empty_list(l2))
		{
		List la,lai,ll;
		lai=NewList();
		la=CommaToList(ListFirst(l2));
		ll=la;
		while(!is_empty_list(ll))
			{
			lai=AppendLast(lai,
				InterfSetAlias(MakeCompound1(OPR_ALIAS,ListFirst(ll)),0));
			ll=ListTail(ll);
			}
		l3=AppendFirst(l3,ProcessAlias(CopyTerm(tt)));
		RemoveAlias(lai);
		RemoveList(la);
		l2=ListTail(l2);
		}
	RemoveList(l1);
	l1=l2plus(l3);
	RemoveList(l3);
	return l1;
	}	



Term WheredTerm(Term t)
	{
	if(is_compound(t) && CompoundName(t)==OPR_WHERE && CompoundArity(t)==2)
		return WheredTerm(namely_w_t(t));
	if(is_atomic(t))
		return t;
	if(is_compound(t))
		{
		Term t1;
		int ac,i;
		ac=CompoundArity(t);
		t1=NewCompound(CompoundFunctor(t));
		for(i=1;i<=ac;i++)
			SetCompoundArg(t1,i,WheredTerm(ConsumeCompoundArg(t,i)));
		FreeAtomic(t);
		return t1;
		}
	if(is_list(t))
		{
		List l,tl;
		tl=t;
		l=NewList();
		while(!is_empty_list(tl))
			{
			l=AppendLast(l,WheredTerm(ListFirst(tl)));
			tl=ListTail(tl);
			}
		RemoveList(t);
		return l;	
		}
	return t;
	}
	
	
	
	
	
	

