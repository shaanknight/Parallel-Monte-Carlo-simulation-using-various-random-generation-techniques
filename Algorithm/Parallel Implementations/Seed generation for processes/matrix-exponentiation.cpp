const int N = 3;
int tmp[N+1],red[N+1][N+1];
void multiply(int pre[N+1][N+1],int base[],ll m,ll mod)
{
	int i,j,k;
	while(m)
	{
		if(m&1ll)
		{
			for(i=1;i<=N;++i)
				tmp[i] = 0;
			for(i=1;i<=N;++i)
				for(j=1;j<=N;++j)
					tmp[i] = (tmp[i] + 1ll*pre[i][j]*base[j]%mod)%mod; 
			for(i=1;i<=N;++i)
				base[i] = tmp[i];
		}
		for(i=1;i<=N;++i)
			for(j=1;j<=N;++j)
				red[i][j] = 0;
		for(i=1;i<=N;++i)
			for(j=1;j<=N;++j)
				for(k=1;k<=N;++k)
					red[i][j] = (red[i][j] + 1ll*pre[i][k]*pre[k][j]%mod)%mod;
		for(i=1;i<=N;++i)
			for(j=1;j<=N;++j)
				pre[i][j] = red[i][j];
		m >>= 1;
	}
}