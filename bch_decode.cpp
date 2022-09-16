#include "mex.h"
#include <iostream>
#include <chrono>
#include <iostream>
#include "math.h"
#include "string"
#include <stdio.h>
#include "time.h"
#include <cstring>

using namespace std;

typedef unsigned char   BYTE;
typedef unsigned long   DWORD;

void Generate_GF(BYTE *PPoly, int PPolyDegree,int m_n,int *m_alpha_to,int *m_index_of)
{
	int	i, mask;
	mask = 1;
	m_alpha_to[PPolyDegree] = 0;
	for(i=0; i<PPolyDegree; i++) 
	{
		m_alpha_to[i] = mask;
		m_index_of[m_alpha_to[i]] = i;
		if (PPoly[i] != 0)
			m_alpha_to[PPolyDegree] ^= mask;
		mask <<= 1;
	}
	m_index_of[m_alpha_to[PPolyDegree]] = PPolyDegree;
	mask >>= 1;
	for(i=PPolyDegree+1; i<m_n; i++) 
	{
		if (m_alpha_to[i-1]>=mask)
			m_alpha_to[i] = m_alpha_to[PPolyDegree] ^ ((m_alpha_to[i - 1] ^ mask) << 1);
		else
			m_alpha_to[i] = m_alpha_to[i-1]<<1;
		m_index_of[m_alpha_to[i]] = i;
	}
	m_index_of[0] = -1;
}

void Field_Init(DWORD FieldGenePol,int *m_alpha_to,int *m_index_of) //域生成多项式的10进制数
{
	int				i;
	BYTE			temp;
	bool			bFlag = false;
	BYTE			PPoly[32], PPolyInv[32];
	int				PPolyDegree;

	memset(PPoly, 0, 32);
	PPolyDegree = 0;
	for(i=31; i>=0; i--)
	{
		temp = (FieldGenePol>>i)&0x01;
		if(temp==1)
		{
			PPoly[PPolyDegree] = 1;
			bFlag = true;
		}
		if(bFlag)
			PPolyDegree++;
	}
	PPolyDegree--;

	int m_n = (1<<PPolyDegree) - 1;

	for(i=0; i<=PPolyDegree; i++)
		PPolyInv[i] = PPoly[PPolyDegree-i];
	Generate_GF(PPolyInv, PPolyDegree,m_n,m_alpha_to,m_index_of);
}


void Decode(BYTE *pCodeBit,DWORD FieldGenePol, int m_t,int  m_n,int m_k, BYTE *pInfoBit)
{


    int *m_alpha_to = new int [m_n];
	memset(m_alpha_to,0,m_n*sizeof(int));

	int *m_index_of = new int [m_n];
	memset(m_index_of,0,m_n*sizeof(int));

	Field_Init(FieldGenePol,m_alpha_to,m_index_of) ;

	int    i, j, u, q, t2, count = 0, syn_error = 0;

	int *d = new int[1026];
	memset(d,0,1026*sizeof(int));

	int *l = new int[1026];
	memset(l,0,1026*sizeof(int));

	int *u_lu = new int[1026];
	memset(u_lu,0,1026*sizeof(int));

	int *s = new int[1025];
	memset(s,0,1025*sizeof(int));

	int *root = new int[200];
	memset(root,0,200*sizeof(int));

	int *loc = new int[200];
	memset(loc,0,200*sizeof(int));

	int *reg = new int[201];
	memset(reg,0,201*sizeof(int));

	int *BCH_elp[1024];

	for (int i = 0;i<1024;i++)
	{
		BCH_elp[i] = new int [1024];
	}




	t2 = 2 * m_t;


	for(i=1; i<=t2; i++)
	{
		s[i] = 0;
		for(j=0; j<m_n; j++)
		{
			if(pCodeBit[j]!=0)
				s[i]^=m_alpha_to[(i*(m_n-1-j))%m_n];
		}
		if(s[i]!=0)
			syn_error = 1; 
		
		s[i] = m_index_of[s[i]];
	}

	if (syn_error) 
	{	
	
		d[0] = 0;			
		d[1] = s[1];		
		BCH_elp[0][0] = 0;		
		BCH_elp[1][0] = 1;		
		for (i = 1; i < t2; i++) 
		{
			BCH_elp[0][i] = -1;	    //存错位多项式的指数形式
			BCH_elp[1][i] = 0;	 //错位多项式
		}
		l[0] = 0;       //初始化D（-1）
		l[1] = 0;       //初始化D（0）
		u_lu[0] = -1;   // j - D(j);
		u_lu[1] = 0;    
		u = 0;

		do
		{
			u++;
			if(d[u]==-1)            // -1 对应 指数 dj = 0 
			{
				l[u + 1] = l[u];    //对应D（j+1） = D(j)
				for(i=0; i<=l[u]; i++) 
				{
					BCH_elp[u + 1][i] = BCH_elp[u][i];        
					BCH_elp[u][i] = m_index_of[BCH_elp[u][i]];
				}
			} 
			else
			{
				q = u - 1;
				while ((d[q] == -1) && (q > 0))   //找i    0< i<j di ！= 0
					q--;      
				
				if(q>0) 
				{
					j = q;
					do 
					{
						j--;
						if ((d[j] != -1) && (u_lu[q] < u_lu[j]))       //求i- D（i) 的最大值
							q = j;
					} while (j > 0);
				}
				
			
				if(l[u]>l[q]+u-q)    //系数  j-i+deerta(i)
					l[u+1] = l[u];
				else
					l[u+1] = l[q]+u-q;   //更新系数 D(j+1)
				
			
				for(i=0; i<t2; i++)
					BCH_elp[u+1][i] = 0;
				for(i=0; i<=l[q]; i++) //最高次为多少，就循环多少次 q = i; u= j； 课本迭代译码中的 i j
				{
					if(BCH_elp[q][i]!=-1)
						BCH_elp[u+1][i+u-q] = m_alpha_to[(d[u]+m_n-d[q]+BCH_elp[q][i])%m_n];  //乘以x的j-i次方，相当于又移动j-i位
				}
				for(i=0; i<=l[u]; i++)  //进行相加 加完后变指数形式
				{
					BCH_elp[u+1][i]^=BCH_elp[u][i];
					BCH_elp[u][i] = m_index_of[BCH_elp[u][i]];
				}
			}
			
			u_lu[u+1] = u + 1 -l[u + 1];  //计算 j - D（j)
			
			
			if(u<t2) 
			{	
				
				if(s[u+1] != -1)
					d[u+1] = m_alpha_to[s[u+1]];     //把指数变10进制数
				else
					d[u+1] = 0;    
				for(i=1; i<=l[u+1]; i++)
				{
					if((s[u+1-i]!=-1)&&(BCH_elp[u+1][i]!=0))
						d[u+1]^= m_alpha_to[(s[u+1-i] + m_index_of[BCH_elp[u+1][i]])%m_n];
				}
					
				d[u+1] = m_index_of[d[u+1]];	
			}
		}while((u<t2)&&(l[u+1]<=m_t));  //j < 2t 且错位多项式的系数小于纠错能力t

		u++;
		if(l[u]<=m_t) 
		{
			for(i=0; i<=l[u]; i++)
				BCH_elp[u][i] = m_index_of[BCH_elp[u][i]];
			
			
			for(i=1; i<=l[u]; i++)
				reg[i] = BCH_elp[u][i];
			count = 0;
			for(i=1; i<=m_n; i++) 
			{
				q = 1;
				int a = 0;
				for(j=1; j<=l[u]; j++)
				{
					if (reg[j]!=-1) 
					{
					

						a = (reg[j]+ j * i)%m_n;
						q^=m_alpha_to[a];

					}
				}
				if (!q) 
				{	
					root[count] = i;
					loc[count] = m_n - i;
					count++;
				}
				
			}
			if (count == l[u])	
				
			{
				for (i = 0; i < l[u]; i++)
					
				    pCodeBit[m_n-1-loc[i]] ^= 1;
			}
			else; 
		}
	}
	memcpy(pInfoBit, pCodeBit, m_k);

	delete []m_alpha_to;
	delete []m_index_of;

	delete []d;
	
	delete []l ;
	
	delete []u_lu;

	delete []s;
	
	delete []root;

	delete []loc ;

    delete [] reg ;
	
	for (int i = 0;i<1024;i++)
	{
		delete []BCH_elp[i];
	}
}


void mexFunction(int nlhs,mxArray *plhs[],int nrhs ,const mxArray *prhs[])  //plhs 输出  prhs 输入
{
	double *input; 
	
	int bchn;
	int bchk;
	int bcht;
	int bchFieldGenePol;

	double *output;

	input = mxGetPr(prhs[0]);
	bchFieldGenePol = *(mxGetPr(prhs[1]));
	bcht = *(mxGetPr(prhs[2]));
	bchn = *(mxGetPr(prhs[3]));
	bchk = *(mxGetPr(prhs[4]));
	

	plhs[0] = mxCreateDoubleMatrix(1,bchk,mxREAL);
	output = mxGetPr(plhs[0]);
	
	
	BYTE *input_byte = (BYTE*)malloc(bchn);
	memset(input_byte,0,bchn);
	
	for(int i = 0 ; i < bchn ;i++)
	{
		if(input[i] == 1)
			input_byte[i] = 1;
	}

	BYTE *output_byte = (BYTE*)malloc(bchk);
	memset(output_byte,0,bchk);

    Decode(input_byte,bchFieldGenePol, bcht,bchn,bchk, output_byte);
  
	for(int i = 0 ; i < bchk ;i++)
	{
		if(output_byte[i] == 1)
			output[i] = 1;
	}

	free(input_byte);
	free(output_byte);

}
