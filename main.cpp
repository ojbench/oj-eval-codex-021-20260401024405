// Optimized O(n log n) solution using segment tree beats with range count
#include <bits/stdc++.h>
using namespace std;

struct SegTree {
    struct Node { int mx, se, cmx; };
    int N; vector<Node> st; 
    int INFV;

    SegTree(int n=0, int infv=0) { if(n) init(n, infv); }
    void init(int n, int infv){ N = n; INFV = infv; st.assign(4*N+4, Node{-1, -1, 0}); build(1,1,N); }
    void build(int p, int l, int r){
        if(l==r){ st[p] = {-1, -1, 1}; return; }
        int m=(l+r)>>1; build(p<<1,l,m); build(p<<1|1,m+1,r); pull(p);
    }
    void pull(int p){
        int lc=p<<1, rc=p<<1|1;
        int mx = max(st[lc].mx, st[rc].mx);
        int se = max(min(st[lc].mx, st[rc].mx), max(st[lc].se, st[rc].se));
        int cmx = 0;
        if(st[lc].mx == mx) cmx += st[lc].cmx;
        if(st[rc].mx == mx) cmx += st[rc].cmx;
        st[p] = {mx, se, cmx};
    }
    void apply_chmin(int p, int x){
        if(st[p].mx <= x) return;
        // Only valid when x > se; typical beats condition should be ensured by caller
        st[p].mx = x;
    }
    void push(int p){
        int lc=p<<1, rc=p<<1|1;
        int lim = st[p].mx;
        if(st[lc].mx > lim) apply_chmin(lc, lim);
        if(st[rc].mx > lim) apply_chmin(rc, lim);
    }
    void range_chmin(int L, int R, int x){ if(L>R) return; range_chmin(1,1,N,L,R,x); }
    void range_chmin(int p, int l, int r, int L, int R, int x){
        if(l>R || r<L || st[p].mx <= x) return;
        if(L<=l && r<=R && st[p].se < x){
            apply_chmin(p, x);
            return;
        }
        push(p);
        int m=(l+r)>>1;
        if(L<=m) range_chmin(p<<1,l,m,L,R,x);
        if(R>m) range_chmin(p<<1|1,m+1,r,L,R,x);
        pull(p);
    }
    void point_set(int pos, int val){ point_set(1,1,N,pos,val); }
    void point_set(int p, int l, int r, int pos, int val){
        if(l==r){
            st[p] = {val, -1, 1};
            return;
        }
        push(p);
        int m=(l+r)>>1;
        if(pos<=m) point_set(p<<1,l,m,pos,val);
        else point_set(p<<1|1,m+1,r,pos,val);
        pull(p);
    }
    int count_ge(int L, int R, int t){ if(L>R) return 0; return count_ge(1,1,N,L,R,t); }
    int count_ge(int p, int l, int r, int L, int R, int t){
        if(l>R || r<L || st[p].mx < t) return 0;
        if(L<=l && r<=R && st[p].se < t) return st[p].cmx;
        push(p);
        int m=(l+r)>>1;
        int res = 0;
        if(L<=m) res += count_ge(p<<1,l,m,L,R,t);
        if(R>m) res += count_ge(p<<1|1,m+1,r,L,R,t);
        return res;
    }
};

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int n; if(!(cin>>n)) return 0;
    vector<pair<long long,long long>> pts(n);
    for(int i=0;i<n;++i) cin>>pts[i].first>>pts[i].second;

    // Compress y to 1..n, sort by x to get permutation
    vector<long long> ys(n);
    for(int i=0;i<n;++i) ys[i]=pts[i].second;
    vector<long long> ys_sorted=ys; sort(ys_sorted.begin(), ys_sorted.end()); ys_sorted.erase(unique(ys_sorted.begin(), ys_sorted.end()), ys_sorted.end());
    vector<int> ord(n); iota(ord.begin(), ord.end(), 0);
    sort(ord.begin(), ord.end(), [&](int a,int b){ if(pts[a].first!=pts[b].first) return pts[a].first<pts[b].first; return pts[a].second<pts[b].second; });
    vector<int> p(n); // permutation of y ranks in x-order
    for(int idx=0; idx<n; ++idx){ int i=ord[idx]; int r = int(lower_bound(ys_sorted.begin(), ys_sorted.end(), pts[i].second) - ys_sorted.begin()) + 1; p[idx]=r; }

    // pos[value] = index position in x-order (1-based)
    vector<int> pos(n+1);
    for(int i=0;i<n;++i) pos[p[i]] = i+1;

    int INF = n+1;
    SegTree st(n, INF);

    long long ans = 0;
    // process values v from 1..n
    for(int v=1; v<=n; ++v){
        if(v>=2){
            int k = v-1;
            int pk = pos[k];
            // activate position pk with NG = INF
            st.point_set(pk, INF);
            // apply chmin on prefix [1, pk-1] with X = pk
            if(pk-1 >= 1) st.range_chmin(1, pk-1, pk);
        }
        int j0 = pos[v];
        // count active positions i<j0 with NG >= j0
        int cnt_ge = st.count_ge(1, j0-1, j0);
        ans += cnt_ge;
    }

    cout << ans << '\n';
    return 0;
}
