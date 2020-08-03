/* Guo, Weilong; guoweilong@gmail.com; 2013-10-19
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <cstring>

using namespace std;

struct parameter
{
	string infile;		// -i
	string outfile;		// -o
	string type;		// -t
	int n;		// -n
	int l;		// -l
	int L;		// -L
	int U;		// -U
};

struct Seq{
	string name;
	string content;
};

parameter param;

void exit_with_help( void )
{
	printf(
		"Usage: GenerateRead_WG_PE -i <input.fa> -n <int> -l <int> -o <output_prefix> [-L 100 -U 300 -t ff/fr/rf/rr]\n"
		"Author: Guo, Weilong; guoweilong@gmail.com; 2013-10-22\n"
		"Description:\n"
		"    Randomly genrate reads from input.fa Whole Genome-widely. Lower/Upper\n"
		"    cases of the original sequence would not be kept in the results.\n"
		"    All reads are from 5\' to 3\' direction.\n"
		"    The chance for one read be generated from one chromosome is due to the\n"
		"    lengths of chromsomes.\n"
		"Options:\n"
		"    -i : input file, Genome sequences file in fasta format\n"
		"    -n : number of read to be generated\n"
		"    -l : length of read to be generated\n"
		"    -o : prefix of the outfile name, in fasta format, all the upper case means methylated\n"
		"    -L : minimum intersection reads, default 100\n"
		"    -U : maximum intersection reads, default 300\n"
		"    -t ff/fr/rf/rr : direction of reads [default: fr]\n"
		"Output format:\n"
		"    >1_chr1_-_611837_611897\n"
		"    cggggctggcaggttcagaggtcgcagcctggagcccagtgagcacagacattagagccg\n"
		"    >2_chr1_+_577257_577317\n"
		"    cggcatcctgcattcctagctgcttcccagagctcccagtgagcacagagccacacctca\n"
		"\n"
	);
	exit(1);
}

vector<string> string_tokenize(const string& str, const string& delimiters = " \t\n\r", bool skip_empty = true);
inline vector<string> string_tokenize(const string& str, const string& delimiters, bool skip_empty) {
	// Skip delimiters at beginning.
	string::size_type lastPos = skip_empty ? str.find_first_not_of(delimiters, 0) : 0;
	// Find first "non-delimiter".
	string::size_type pos     = str.find_first_of(delimiters, lastPos);
	vector<string> result;
	result.clear();

	while (string::npos != pos || string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		//__ASSERT(pos > lastPos || !skip_empty, "internal error, pos <= lastPos.\n");

		//if (pos == lastPos) result.push_back("");
		result.push_back(str.substr(lastPos, pos - lastPos));

		if (pos == string::npos) break;
		if (pos == str.length() - 1) {
			if (!skip_empty) result.push_back("");
			break;
		}
		// Skip delimiters.  Note the "not_of"
		lastPos = skip_empty ? str.find_first_not_of(delimiters, pos) : pos + 1;
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
	return result;
}

void ToUpperString(string &str)
{
	transform(str.begin(), str.end(), str.begin(), (int (*)(int))toupper);
}

void ToLowerString(string &str)
{
	transform(str.begin(), str.end(), str.begin(), (int (*)(int))tolower);
}

void parse_command_line(int argc, char **argv)
{
	// GenomePosition -f <filename.fa> -s TACAG -A -o <output.fa>
	int i;

	for(i=2;i<argc;i++)
	{
		if(argv[i-1][0] != '-') break;
		switch(argv[i-1][1])
		{
			case 'i':
				param.infile = string(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'o':
				param.outfile = string(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 't':
				param.type = string(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'n':
				param.n = atoi(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'l':
				param.l = atoi(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'L':
				param.L = atoi(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			case 'U':
				param.U = atoi(argv[i]);
				if(++i>argc)	exit_with_help();
				break;
			default:
				fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
		}
	}
	if ( !param.infile.length() || !param.outfile.length() ){
		exit_with_help();
	}
}

void init(){
	// default values
	param.infile = "";
	param.outfile = "";
	param.n = 0;
	param.l = 20;
	param.L = 100;
	param.U = 300;
}

// Read the fasta file of genome, storing in [Genome] structure
int ReadGenome( map<string, string> & Genome )
{
	ifstream FAfile( param.infile.c_str() );
	if(!FAfile) {
		printf("cannot open geneNome file %s\n", param.infile.c_str() );
		return -1;
	}
	string chr = "";
	string chr_long_strings = "";
	while(!FAfile.eof()){
		char buffer[10000 + 1];
		FAfile.getline(buffer, 10000);
		if(buffer[strlen(buffer) - 1] == '\r'){
			buffer[strlen(buffer) - 1] = '\0';
		}
		string tmp = buffer;
		if(buffer[0] == '>'){
			if(chr == ""){
				chr_long_strings = "";
				vector<string> tokens = string_tokenize(tmp.substr(1));
				chr = tokens[0];	//specify which chromosome it is
				continue;
			}
			Genome[chr] = chr_long_strings;
			chr_long_strings = "";
			vector<string> tokens = string_tokenize(tmp.substr(1));
			chr = tokens[0];
		}else{
			chr_long_strings += tmp;
		}
	}
	Genome[chr] = chr_long_strings;
	// store the sequence in chr2seq[chr]
	FAfile.close();
	return 1;
}

// Output string in lines, length of each line <= "step"
int OutputString ( string & str, int step ) {
	int len = str.size();
	for (int i=0; i < len -1; i+=step) {
		cout << str.substr(i, step) << endl;
	}
	return 0;
}

// Output the genome as fasta format
int OutputGenome( map<string, string> & Genome )
{
	map<string, string>::iterator giter;
	for ( giter=Genome.begin(); giter!=Genome.end(); giter++ ) {
		string chr = giter->first;
		string chr_long_str = Genome[chr];
		cout << ">" << chr << endl;
		OutputString (chr_long_str, 50);
	}
	return 0;
}

// Get the reversed complementary sequence, 5'->3'
string Antisense ( string read ) {
	int len = read.length();
	for (int i = 0; i < len; i++ ) {
		switch(read[i]){
		case 'a': read[i]='t'; break;
		case 'A': read[i]='T'; break;
		case 'c': read[i]='g'; break;
		case 'C': read[i]='G'; break;
		case 'g': read[i]='c'; break;
		case 'G': read[i]='C'; break;
		case 't': read[i]='a'; break;
		case 'T': read[i]='A'; break;
		default: break;
		}
	}
	// Get the reverse of string
	return  string ( read.rbegin(), read.rend() );
}

// Find the right chromosome and position on that chromosome by random int
int FindPosition( int max, vector<int> & chrcummlen, vector<string> & chrkey, string & randchr,
				  int & randpos, map<string, string> & Genome ) 
{
	int n = rand() % max;
	for (int i=0; i<chrcummlen.size(); i++) {
		if ( n < chrcummlen[i] ) {
			randchr = chrkey[i];
			int len = Genome[randchr].length();
			randpos = rand() % len;
			return 1; // Success
		}
	}
	return 0; // Fail
}

int GetType(string type) {
	if ( type==string("ff") )
		return 0;
	if ( type==string("fr") )
		return 1;
	if ( type==string("rf") )
		return 2;
	if ( type==string("rr") )
		return 3;
	return 1;
}

// Generate Reads randomly
int GenerateReads ( map<string, string> & Genome ) 
{
	ofstream ofile_1( (param.outfile + "_1.fa").c_str() );
	ofstream ofile_2( (param.outfile + "_2.fa").c_str() );
	srand(99); // Important, set a rand seed here
	// Store all keys into a vector
	vector<string> chrkey;
	vector<int> chrcummlen; // chromosome cummulated length
	int cummlen = 0;
	cerr << "If output numbers contain negative one, your genome size is ";
	cerr << "too large and make program overfloated" << endl;
	for ( map<string, string>::iterator giter = Genome.begin(); giter != Genome.end(); giter++ ) {
		chrkey.push_back(giter->first);
		cerr << giter->first << endl;
		cummlen += giter->second.length()/100; // ex. mouse genome size will get INT overflow
		chrcummlen.push_back(cummlen);
		cerr << cummlen << endl;
	}
	int type = GetType(param.type);
	//cout << "type=" << type << endl;
	int nchr = chrkey.size();
	// Randomly generate the reads
	for ( int i = 0; i < param.n; i++) {
		//cout << "NO:" << i << endl;
		string randchr;
		int randpos;
		string read, read_1, read_2;
		int dir;
		int len;
		bool state = true;
		while (state)
		{
			while ( !FindPosition( cummlen, chrcummlen, chrkey, randchr, randpos, Genome ) )
				;
			dir = rand()%2; // direction
			double p = (double)rand() / RAND_MAX;
			len = int(param.l * 2 + param.L * p + param.U * (1-p));
			if ( randpos + len < Genome[randchr].length() ) {
				read = Genome[randchr].substr(randpos, len);
				if (dir == 1)
					read = Antisense(read);
				switch (type) {
					case 0 : // ff
						read_1 = read.substr(0, param.l);
						read_2 = read.substr(len - param.l, param.l);
						break;
					case 1: // fr
						read_1 = read.substr(0, param.l);
						read_2 = Antisense(read).substr(0, param.l);
						break;
					case 2: // rf
						read_1 = Antisense(read.substr(0, param.l));
						read_2 = read.substr(len - param.l, param.l);
						break;
					case 3: // rr
						read_1 = Antisense(read.substr(0, param.l));
						read_2 = Antisense(read).substr(0, param.l);
						break;
					default:
						break;
				}
				//cout << read_1 << endl;
				//cout << read_2 << endl;
				if ( read.find('n')==string::npos && read.find('N')==string::npos )
				{
					state = false; // an out gate
				}
			}
		}
		//cout << dir << endl;
		if ( dir == 0 ) { // plus strand
			// First end
			ofile_1 << ">" << i << "_" << randchr;
			if ( type == 0 || type == 1) {
				ofile_1 << "_+_";
			}else{
				ofile_1 << "_-_";
			}
			ofile_1 << randpos << "_" << randpos + param.l << endl;
			// Second end
			ofile_2 << ">" << i << "_" << randchr; 
			if ( type == 0 || type == 2) {
				ofile_2 << "_+_"; 
			} else {
				ofile_2 << "_-_";
			} 
			ofile_2 << randpos + len - param.l << "_" << randpos + len << endl;
		} else { // minus strand
			// First end
			ofile_1 << ">" << i << "_" << randchr;
			if ( type == 0 || type == 1) {
				ofile_1 << "_-_";
			} else {
				ofile_1 << "_+_";
			} 
			ofile_1 << randpos + len - param.l << "_" << randpos + len << endl;
			// Second end
			ofile_2 << ">" << i << "_" << randchr;
			if ( type == 0 || type == 2) {
				ofile_2 << "_-_"; 
			} else {
				ofile_2 << "_+_";
			} 
			ofile_2 << randpos << "_" << randpos + param.l << endl;
		}
		ofile_1 << read_1 << endl;
		ofile_2 << read_2 << endl;
	}
	ofile_1.close();
	ofile_2.close();
	return 0;
}


int main(int argc, char* argv[])
{
	init();
	parse_command_line(argc,argv);

	map<string, string> Genome;
	ReadGenome ( Genome );
//	OutputGenome ( Genome );
	GenerateReads ( Genome );
	return 1;
}

