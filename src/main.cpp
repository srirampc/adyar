//#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <sdsl/suffix_arrays.hpp>
#include <sdsl/lcp_bitcompressed.hpp>
#include <sdsl/rmq_support.hpp>
#include <time.h>
#include <vector>

#include "util.hpp"

struct RunArgs{
	std::string x;
	std::string y;
	int k;
};

struct GST{
	sdsl::csa_bitcompressed<> SA;
	sdsl::lcp_bitcompressed<> LCP;
	std::vector<int> ISA;
	std::string text;
};

void read_input(int argc, char *argv[],
				RunArgs& args) {
	std::string fname(argv[1]);
	args.k = atoi(argv[2]);

	std::ifstream file_object;
	file_object.open(fname);
	std::string text1 = "";
	std::string text2;
	file_object >> text1;
	file_object >> text2;
	file_object.close();

	args.x = trim(text1);
	args.y = trim(text2);

	std::cout << "# mismatches (k)  : " << args.k << std::endl;
	std::cout << "X length          : " << args.x.length() << std::endl;
	std::cout << "Y length          : " << args.y.length() << std::endl;
}

void construct_gst(std::string& sx, std::string& sy, GST& gst){
	gst.text = sx + std::string("$") + sy; // + std::string("@");
	construct_im(gst.SA, gst.text, 1);
	construct_im(gst.LCP, gst.text, 1);
	gst.ISA.resize(gst.SA.size());

	for(auto i=0u; i < gst.SA.size(); i++) {
		gst.ISA[gst.SA[i]] = i;
	}

	std::cout << "gst.text.length()  : " << gst.text.length() << std::endl;
	std::cout << "gst.SA.size()      : " << gst.SA.size() << std::endl;
	std::cout << "gst.LCP.size()     : " << gst.LCP.size() << std::endl;
	std::cout << "gst.ISA.size()     : " << gst.ISA.size() << std::endl;
}

void find_matches(RunArgs& args, GST& gst, 
                  std::vector<int>& match_length,
				  std::vector<int>& left_match,
				  std::vector<int>& right_match){
    int n = gst.text.length() + 1;
	int m = args.y.length();
	match_length.resize(n, 0);
	left_match.resize(n, 0);
	right_match.resize(n, 0);
	int min=0, last_match=0;
	for(int i=1; i < n-1; i++) {
		if ((gst.SA[i] > n-m-2 && gst.SA[i+1] > n-m-2) || (gst.SA[i] < n-m-2 && gst.SA[i+1] < n-m-2)) {

			if (gst.LCP[i+1] <= min)
				min = gst.LCP[i+1];
			match_length[i+1] = min;
			left_match[i+1] = last_match;
		}
		else {

			min = gst.LCP[i+1];
			match_length[i+1] = gst.LCP[i+1];
			left_match[i+1] = i; //stores the position in second string where the match occures
			last_match = i;
		}
	}

	min=0; last_match=0;
	for(int i=n-1; i >= 1; i--) {
		if ((gst.SA[i] > n-m-2 && gst.SA[i-1] > n-m-2) || (gst.SA[i] < n-m-2 && gst.SA[i-1] < n-m-2)) {

			if (gst.LCP[i] <= min)
				min = gst.LCP[i];
			if(min>match_length[i-1]){
				match_length[i-1] = min;
				right_match[i-1] = last_match;
				left_match[i-1] = 0;
			}
			else if(min<match_length[i-1])
				right_match[i-1]=0;
			else
				right_match[i-1]=last_match;

		}
		else {

			min = gst.LCP[i];
			last_match = i;
			if(min>match_length[i-1]){
				match_length[i-1] = min;
				right_match[i-1] = last_match;
				left_match[i-1] = 0;
			}
			else if(min<match_length[i-1])
				right_match[i-1] = 0;
			else
				right_match[i-1] = last_match;
		}
	}
}

void compute_acs_zero(RunArgs& args, GST& gst, 
              std::vector<int>& match_length){
    int n = gst.text.length() + 1;
	int m = args.y.length();
	std::vector<int> s1(n-m-2, 0);
	std::vector<int> s2(m, 0);

	// break s-array into individual s-i-array for each text-i and
	// reorder to lexographical order 
	for(int i=2; i<n; i++) { // start from 2 since first two are bad.
		if (gst.SA[i] >= n-m-1)
			s2[gst.SA[i]-(n-m-1)] = match_length[i];
		if (gst.SA[i] <= n-m-3)
			s1[gst.SA[i]] = match_length[i];
	}
	// calculate avgs1 and avgs2 and d_acs
	double avg_s1=0, avg_s2=0;
	for(int i=0; i<n-m-2; i++) {avg_s1 += s1[i];}
	avg_s1 = avg_s1/(n-m-2);
	for(int i=0; i<m; i++) {avg_s2 += s2[i];}
	avg_s2 = avg_s2/m;

	double d_acs=0;
	d_acs = ((log10(n-m-2)/(2*avg_s2)) + (log10(m)/(2*avg_s1))) - ((log10(n-m-2)/((n-m-2))) + (log10(m)/(m)));
	std::cout << "ACS for k=0 : " << d_acs << std::endl;
}


int forward_match(std::string const& T, int pos1, int pos2, int k, int m, int n){
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

void compute_acsk(RunArgs& args, GST& gst, 
				  std::vector<int>& match_length,
				  std::vector<int>& left_match,
				  std::vector<int>& right_match){
    int n = gst.text.length() + 1;
	int m = args.y.length();
	std::vector<int> s1(n-m-2, 0);
	std::vector<int> s2(m, 0);
	std::vector<int> sk(n, 0);
	int pos1, p, max=0, pos2=0, flag=0, lcp;

	// calculate sk-array
	for(int i=0; i<n; i++) {
		sk[i] = match_length[i];
		pos1 = gst.SA[i]+match_length[i]+1;
		pos2 = 0;
		max = 0;
		flag=0;

		// forward matching char by char
		if(gst.SA[i] < n-m-2){
			if(left_match[i]>0 && pos1 < n-m-2){
				p = left_match[i];
				while(gst.LCP[p+1] >= match_length[i] && p>0){
					if(gst.SA[p] > n-m-2){
						lcp = forward_match(gst.text, pos1-1, gst.SA[p]+match_length[i], args.k, n-m-2, n);
						if(lcp>max){
							max=lcp;
							pos2 = p;
							flag = 1;
						}
					}
					p--;
				}
			}
			if(right_match[i]>0 && pos1 < n-m-2){
				p = right_match[i];
				while(gst.LCP[p] >= match_length[i] && p<n){
					if(gst.SA[p] > n-m-2){
						lcp = forward_match(gst.text, pos1-1, gst.SA[p]+match_length[i], args.k, n-m-2, n);
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
				left_match[i] = pos2;
				right_match[i] = 0;
			}
			else if(flag == 2){
				left_match[i] = 0;
				right_match[i] = pos2;
			}
			sk[i] = match_length[i]+max;
		}
		else if(gst.SA[i] > n-m-2){
			if(left_match[i]>0 && pos1 < n){
				p = left_match[i];
				while(gst.LCP[p+1] >= match_length[i] && p>0){
					if(gst.SA[p] < n-m-2){
						lcp = forward_match(gst.text, gst.SA[p]+match_length[i], pos1-1, args.k, n-m-2, n);
						if(lcp>max){
							max=lcp;
							pos2 = p;
							flag = 1;
						}
					}
					p--;
				}
			}
			if(right_match[i]>0 && pos1 < n){
				p = right_match[i];
				while(gst.LCP[p] >= match_length[i] && p<n){
					if(gst.SA[p] < n-m-2){
						lcp = forward_match(gst.text, gst.SA[p]+match_length[i], pos1-1, args.k, n-m-2, n);
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
				left_match[i] = pos2;
				right_match[i] = 0;
			}
			else if(flag == 2){
				left_match[i] = 0;
				right_match[i] = pos2;
			}
			sk[i] = match_length[i]+max;
		}

	}

	for(int i=1; i<n; i++) {
		if(left_match[i] < right_match[i])
			left_match[i] = right_match[i];
		//cout << i << " "<< gst.SA[left_match[i]] << endl;
	}

	// backward and forward search and storage of missmatches for each i
	int missmatch_position[2*args.k+2];
	for(int i=0;i<(2*args.k)+2;i++){missmatch_position[i] = -1;}
	int k_counter;
	int j=0, substring_len=0;
	//cout << "doing matching..." << endl;
	for(int i=1; i<n; i++) { // skip i = 0
		// initialize pos1 , pos2
		pos1 = gst.SA[i]; pos2 = gst.SA[left_match[i]];
		k_counter = args.k;
		//backward matching
		//cout << "entering backward while loop at " << i << endl;
		while(((pos1 < n-m-2 && pos2 > n-m-2) || (pos1 > n-m-2 && pos2 < n-m-2)) && k_counter >= 0 && pos1 >= 0 && pos2 >= 0) {
			if(gst.text[pos1] == gst.text[pos2]) {
				pos1--; pos2--;
			}
			else {
				missmatch_position[k_counter] = pos1+1; pos1--; pos2--; k_counter--;
			}
		}
		//forward matching
		pos1 = gst.SA[i] + match_length[i] + 1; pos2 = gst.SA[left_match[i]] + match_length[i] + 1;
		k_counter = args.k;
		//cout << "entering forward while loop at " << i << endl;
		while(((pos1 < n-m-2 && pos2 > n-m-2) || (pos1 > n-m-2 && pos2 < n-m-2)) && k_counter >= 0 && pos1 < n && pos2 < n) {
			if(gst.text[pos1] == gst.text[pos2]) {
				pos1++; pos2++;
			}
			else {
				missmatch_position[2*args.k -k_counter+1] = pos1-1; k_counter--;
			}
		}
		// for each missmatch compair sk[i] and [k+j+1]th - [j]th
		// elements of missmatch_position array replace if less
		j=0;
		for(j=0;j<=args.k;j++){
			if(missmatch_position[j] == -1 || missmatch_position[j+args.k+1] == -1)
				continue;
			substring_len = missmatch_position[args.k+j+1] - missmatch_position[j]+1;
			if(substring_len>sk[missmatch_position[j]])
				sk[missmatch_position[j]] = substring_len;
		}
		j=0;
		for(j=0;j<(2*args.k)+2;j++){missmatch_position[j] = -1;}
	}
	//cout << "done with heavy stuff" << endl;
	// modify sk-array so that sk[i] = x implies sk[i+1] >= x-1
	for(int i=1;i<n;i++){
		if(sk[i] > sk[gst.ISA[gst.SA[i] + 1]])
			sk[gst.ISA[gst.SA[i] + 1]] = sk[i] - 1;
	}

	for(int i=1; i<n-1; i++) {
		if (gst.SA[i] >= n-m-1)
			s2[gst.SA[i]-(n-m-1)] = sk[i];
		if (gst.SA[i] <= n-m-3)
			s1[gst.SA[i]] = sk[i];
	}
	//for(i=0; i<n-m-2; i++) {cout << s1[i] << endl;}
	//for(i=0; i<m; i++) {cout << s2[i] << endl;}
	// calculate avgs1 and svgs2 and d_acs
	double avg_s1=0, avg_s2=0;
	for(int i=0; i<n-m-2; i++) {avg_s1 += s1[i];}
	avg_s1 = avg_s1/(n-m-2);
	for(int i=0; i<m; i++) {avg_s2 += s2[i];}
	avg_s2 = avg_s2/m;

	double d_acs=0;
	d_acs = ((log10(n-m-2)/(2*avg_s2)) + (log10(m)/(2*avg_s1))) - ((log10(n-m-2)/((n-m-2))) + (log10(m)/(m)));
	std::cout << "ACS for k= " << args.k << " : " << d_acs << std::endl;
}

int main(int argc, char *argv[]) {
	//cout << text1 << endl;
	//cout << text2 << endl;
	RunArgs args;
	read_input(argc, argv, args);

    GST gst;
	construct_gst(args.x, args.y, gst);

	time_t t;
	t = clock();

	std::vector<int> match_length, left_match, right_match;

    // find the longest matching suffixes to the left and right 
	//  of suffix in GST
    find_matches(args, gst, match_length, left_match, right_match);

	if(args.k==0){ // ACS for k = 0 can be computed from the match lengths
        compute_acs_zero(args, gst, match_length);
	}
	else {
        compute_acsk(args, gst, match_length, left_match, right_match);
	}
	
	t = clock() - t;
	std::cout << "time = " << (float)t/CLOCKS_PER_SEC << " (s) " << std::endl;
	return 0;
}
