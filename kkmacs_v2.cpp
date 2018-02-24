//#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/lcp_bitcompressed.hpp>
#include <sdsl/rmq_support.hpp>
#include <time.h>

using namespace sdsl;
using namespace std;

int forward_match(string const& T, int pos1, int pos2, int k, int m, int n){
	int lcp = 0;
	while(k>=0 && pos1< n-m-2 && pos2 < n){
		if(T[pos1] == T[pos2]){
			lcp++; pos1++; pos2++;
		}
		else {
			pos1++; pos2++; k--;
		}
	}
	return lcp;
}

int main(int argc, char *argv[]) {
  ifstream file_object;
  file_object.open(argv[1]);
  string text1 = "";
  string text2;
  file_object >> text1;
  file_object >> text2;

  //cout << text1 << endl;
  //cout << text2 << endl;

  file_object.close();
  int n = text1.length() + 1;
	int m = text2.length();
	csa_bitcompressed<> SA;
  lcp_bitcompressed<> LCP;
  construct_im(SA, text1, 1);
  construct_im(LCP, text1, 1);
  cout << "SA.size(): " << SA.size() << endl;

	int s[n]={0};
	int min = 0;
	int i;
	int SA_inverse[n];

	for(i=0; i<n; i++){
		SA_inverse[SA[i]] = i;
	}

	int s_position[n]={0};
	int s_position2[n]={0};
	int last_match=0;
	time_t t;
	t = clock();
	for(i=1; i < n-1; i++) {
		if ((SA[i] > n-m-2 && SA[i+1] > n-m-2) || (SA[i] < n-m-2 && SA[i+1] < n-m-2)) {

			if (LCP[i+1] <= min)
				min = LCP[i+1];
			s[i+1] = min;
			s_position[i+1] = last_match;
		}
		else {

			min = LCP[i+1];
			s[i+1] = LCP[i+1];
			s_position[i+1] = i; //stores the position in second string where the match occures
			last_match = i;
		}
	}

	min=0; last_match=0;
	for(i=n-1; i >= 1; i--) {
		if ((SA[i] > n-m-2 && SA[i-1] > n-m-2) || (SA[i] < n-m-2 && SA[i-1] < n-m-2)) {

			if (LCP[i] <= min)
				min = LCP[i];
			if(min>s[i-1]){
				s[i-1] = min;
				s_position2[i-1] = last_match;
				s_position[i-1] = 0;
			}
			else if(min<s[i-1])
				s_position2[i-1]=0;
			else
				s_position2[i-1]=last_match;

		}
		else {

			min = LCP[i];
			last_match = i;
			if(min>s[i-1]){
				s[i-1] = min;
				s_position2[i-1] = last_match;
				s_position[i-1] = 0;
			}
			else if(min<s[i-1])
				s_position2[i-1] = 0;
			else
				s_position2[i-1] = last_match;
		}
	}

	int s1[n-m-2] = {0};
	int s2[m] = {0};

	int k = atoi(argv[2]);

	if(k==0){

	/* break s-array into individual s-i-array for each text-i and
	reorder to lexographical order */
		for(i=0; i<n; i++) {
			if (SA[i] >= n-m-1)
				s2[SA[i]-(n-m-1)] = s[i];
			if (SA[i] <= n-m-3)
				s1[SA[i]] = s[i];
		}
		// calculate avgs1 and avgs2 and d_acs
		double avg_s1=0, avg_s2=0;
		for(i=0; i<n-m-2; i++) {avg_s1 += s1[i];}
		avg_s1 = avg_s1/(n-m-2);
		for(i=0; i<m; i++) {avg_s2 += s2[i];}
		avg_s2 = avg_s2/m;

		double d_acs=0;
		d_acs = ((log10(n-m-2)/(2*avg_s2)) + (log10(m)/(2*avg_s1))) - ((log10(n-m-2)/((n-m-2))) + (log10(m)/(m)));
		printf("ACS for k=0: %f\n", d_acs);
	}

	else {
		int sk[n] = {0};
		int pos1, p, max=0, pos2=0, flag=0, lcp;

			// calculate sk-array
		for(i=0; i<n; i++) {
			sk[i] = s[i];
			pos1 = SA[i]+s[i]+1;
			pos2 = 0;
			max = 0;
			flag=0;

				// forward matching char by char
			if(SA[i] < n-m-2){
				if(s_position[i]>0 && pos1 < n-m-2){
					p = s_position[i];
					while(LCP[p+1] >= s[i] && p>0){
						if(SA[p] > n-m-2){
							lcp = forward_match(text1, pos1-1, SA[p]+s[i], k, n-m-2, n);
							if(lcp>max){
								max=lcp;
								pos2 = p;
								flag = 1;
							}
						}
						p--;
					}
				}
				if(s_position2[i]>0 && pos1 < n-m-2){
					p = s_position2[i];
					while(LCP[p] >= s[i] && p<n){
						if(SA[p] > n-m-2){
							lcp = forward_match(text1, pos1-1, SA[p]+s[i], k, n-m-2, n);
							if(lcp>max){
								max=lcp;
								pos2 = p;
								flag=2;
							}
						}
						p++;
					}
				}
				if(flag == 1) {
					s_position[i] = pos2;
					s_position2[i] = 0;
				}
				else if(flag == 2){
					s_position[i] = 0;
					s_position2[i] = pos2;
				}
				sk[i] = s[i]+max;
			}
			else if(SA[i] > n-m-2){
				if(s_position[i]>0 && pos1 < n){
					p = s_position[i];
					while(LCP[p+1] >= s[i] && p>0){
						if(SA[p] < n-m-2){
							lcp = forward_match(text1, SA[p]+s[i], pos1-1, k, n-m-2, n);
							if(lcp>max){
								max=lcp;
								pos2 = p;
								flag = 1;
							}
						}
						p--;
					}
				}
				if(s_position2[i]>0 && pos1 < n){
					p = s_position2[i];
					while(LCP[p] >= s[i] && p<n){
						if(SA[p] < n-m-2){
							lcp = forward_match(text1, SA[p]+s[i], pos1-1, k, n-m-2, n);
							if(lcp>max){
								max=lcp;
								pos2 = p;
								flag=2;
							}
						}
						p++;
					}
				}
				if(flag == 1) {
					s_position[i] = pos2;
					s_position2[i] = 0;
				}
				else if(flag == 2){
					s_position[i] = 0;
					s_position2[i] = pos2;
				}
				sk[i] = s[i]+max;
			}

			}

		for(i=1; i<n; i++) {
			if(s_position[i] < s_position2[i])
				s_position[i] = s_position2[i];
			//cout << i << " "<< SA[s_position[i]] << endl;
		}

		// backward and forward search and storage of missmatches for each i
		int missmatch_position[2*k+2];
		for(i=0;i<(2*k)+2;i++){missmatch_position[i] = -1;}
		int k_counter;
		int j=0, substring_len=0;
		//cout << "doing matching..." << endl;
		for(i=1; i<n; i++) { // skip i = 0
			// initialize pos1 , pos2
			pos1 = SA[i]; pos2 = SA[s_position[i]];
			k_counter = k;
			//backward matching
			//cout << "entering backward while loop at " << i << endl;
			while(((pos1 < n-m-2 && pos2 > n-m-2) || (pos1 > n-m-2 && pos2 < n-m-2)) && k_counter >= 0 && pos1 >= 0 && pos2 >= 0) {
				if(text1[pos1] == text1[pos2]) {
					pos1--; pos2--;
				}
				else {
					missmatch_position[k_counter] = pos1+1; pos1--; pos2--; k_counter--;
				}
    	}
			//forward matching
			pos1 = SA[i] + s[i] + 1; pos2 = SA[s_position[i]] + s[i] + 1;
			k_counter = k;
			//cout << "entering forward while loop at " << i << endl;
			while(((pos1 < n-m-2 && pos2 > n-m-2) || (pos1 > n-m-2 && pos2 < n-m-2)) && k_counter >= 0 && pos1 < n && pos2 < n) {
				if(text1[pos1] == text1[pos2]) {
					pos1++; pos2++;
				}
				else {
					missmatch_position[2*k -k_counter+1] = pos1-1; k_counter--;
				}
    	}
			/* for each missmatch compair sk[i] and [k+j+1]th - [j]th
			elements of missmatch_position array replace if less*/
			j=0;
			for(j=0;j<=k;j++){
				if(missmatch_position[j] == -1 || missmatch_position[j+k+1] == -1)
					continue;
				substring_len = missmatch_position[k+j+1] - missmatch_position[j]+1;
				if(substring_len>sk[missmatch_position[j]])
					sk[missmatch_position[j]] = substring_len;
			}
			j=0;
			for(j=0;j<(2*k)+2;j++){missmatch_position[j] = -1;}
		}
		//cout << "done with heavy stuff" << endl;
		// modify sk-array so that sk[i] = x implies sk[i+1] >= x-1
		for(i=1;i<n;i++){
			if(sk[i] > sk[SA_inverse[SA[i] + 1]])
				sk[SA_inverse[SA[i] + 1]] = sk[i] - 1;
		}

		for(i=1; i<n-1; i++) {
			if (SA[i] >= n-m-1)
				s2[SA[i]-(n-m-1)] = sk[i];
			if (SA[i] <= n-m-3)
				s1[SA[i]] = sk[i];
		}
		//for(i=0; i<n-m-2; i++) {cout << s1[i] << endl;}
		//for(i=0; i<m; i++) {cout << s2[i] << endl;}
		// calculate avgs1 and svgs2 and d_acs
		double avg_s1=0, avg_s2=0;
		for(i=0; i<n-m-2; i++) {avg_s1 += s1[i];}
		avg_s1 = avg_s1/(n-m-2);
		for(i=0; i<m; i++) {avg_s2 += s2[i];}
		avg_s2 = avg_s2/m;

		double d_acs=0;
		d_acs = ((log10(n-m-2)/(2*avg_s2)) + (log10(m)/(2*avg_s1))) - ((log10(n-m-2)/((n-m-2))) + (log10(m)/(m)));
		t = clock() - t;
		cout << "time = " << (float)t/CLOCKS_PER_SEC << endl;
		printf("ACS for k=%d: %f\n",atoi(argv[2]), d_acs);

	}
  return 0;

}
