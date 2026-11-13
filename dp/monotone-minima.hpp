#pragma once
// f(row, old_column, new_column): true iff the new column is no worse.
// Row argmins must be nondecreasing. Rightmost argmin on ties.
template <typename F> std::vector<int> monotone_minima(int h,int w,const F& f) {
    std::vector<int> dp(h,-1);
    if(h==0||w==0) return dp;
    auto rec=[&](auto&& self,int u,int d,int l,int r)->void {
        if(u>=d) return;
        int m=u+(d-u)/2,best=l;
        for(int j=l+1;j<r;++j) if(f(m,best,j)) best=j;
        dp[m]=best;
        self(self,u,m,l,best+1);
        self(self,m+1,d,best,r);
    };
    rec(rec,0,h,0,w); return dp;
}
