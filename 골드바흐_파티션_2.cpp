#define _USE_MATH_DEFINES
#include <iostream>
#include <complex>
#include <algorithm>
#include <vector>
using namespace std;

typedef complex<double> cpx;
const int SZ = 1 << 21;

void FFT(cpx f[], bool inv = false){
    int n = SZ;
    for(int i = 1, j = 0; i < n; ++i){
        int b = n/2;
        while(!((j ^= b) & b)) b /= 2;
        if(i < j) swap(f[i], f[j]);
    }
    for(int k = 1; k < n; k *= 2){
        double a = (inv ? M_PI/k : -M_PI/k);
        cpx w(cos(a), sin(a));
        for(int i = 0; i < n; i += k*2){
            cpx wp(1, 0);
            for(int j = 0; j < k; ++j){
                cpx x = f[i+j], y = f[i+j+k] * wp;
                f[i+j] = x + y;
                f[i+j+k] = x - y;
                wp *= w;
            }
        }
    }
    if(inv){
        for(int i = 0; i < n; ++i)
            f[i] /= n;
    }
}

void pw(cpx f[]){
	FFT(f);
	for(int i=0; i<SZ; i++) f[i] *= f[i];
	FFT(f, 1);
}

int chk[SZ >> 1] = {1 ,1};
cpx f[SZ];

int main(){
	ios_base::sync_with_stdio(0); cin.tie(0);
	chk[0] = chk[1] = 1;
	for(long long i=2; i<1000000; i++){
		if(chk[i]) continue;
		f[i] = cpx(1, 0);

		for(long long j=i*i; j<1000000; j+=i) chk[j] = 1;
	}
	pw(f);

	int t; cin >> t;
	while(t--){
		int x; cin >> x;
		int ans = round(f[x].real());
		int flag = !chk[x/2];
		if(x/2%2 == 0 && x/2 > 2) flag = 0;
		if(flag) ans--;
		cout << ans/2+flag << "\n";
	}
}