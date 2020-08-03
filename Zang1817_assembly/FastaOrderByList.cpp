/* Usage: FastaOrderByList [-i <input.fa>]  -o <output.fa> -c <ChrList>
 * Weilong GUO, start from 2018-03-15
 * *********************
 * Modification log:
 */

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>

using namespace std;
#include "math.h"

struct parameter
{
	string infile;		// -i
	string outfile;		// -o
	string chrlist;		// -c
	int    len;		// -l
	bool   Info;		// -I
};

struct ChrDir
{
	string chr;
	char dir;
};


class Seq{
public:
	string name;
	string content;
	Seq( void ) {
		name = "";
		content = "";
	}
	Seq(string NAME, string CONTENT) {
		name = NAME;
		content = CONTENT;
	}
};

parameter param;

void exit_with_help( void )
{
	printf(
		"Usage: FastaOrderByList [-i <input.fa>]  -o <output.fa> -c <ChrList>\n"
		"                   [-l 60 -I]\n"
		"Author: Guo, Weilong; guoweilong@126.com; 2018-03-15\n"
		"Last updated: 2018-04-02\n"
		"Description: ReOrdering fasta entries by ChrList\n"
		"Options:\n"
		"  -i <STRING> : Input file, Fasta, STDIN if omitted\n"
		"  -o <STRING> : Output file, Fasta, constructed new file\n"
		"  -c <STRING> : Input Chr list file, TXT, one chr one line, in order of purpose\n"
		"  -l <INT> : Max length of line for output sequence [Default: 60]\n"
		"  -I : Output interactive information when specified\n"
		"  -h : help\n"
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


map<string, string> chr2seq;

/** Output the structure in a FASTA-format file. By default, the length of each
 *  line is 50 bp.
 */
int OutputSequence( ofstream & of, string & seq, int step = 50 )
{
	if(!of) {
		cerr << "Error with the fasta output handle.\n";
		return -1;
	}
	int read_len = seq.length();
	int max_len = read_len - step;
	int i;
	for (i=0; i<max_len; i+=step) {
		of << seq.substr(i, step) << endl;
	}
	if ( i<=(read_len-1) ) {
		of << seq.substr(i, read_len-i) << endl;
	}
	return seq.length();
}



int ReadChrList ( string infile, vector<ChrDir> & ChrList ) {
	ifstream IN(infile.c_str());
	if(!IN) {
		cout << "cannot open input file" << infile.c_str() << endl;
		return -1;
	}
	while ( !IN.eof() ) {
		string str;
		getline(IN, str);
		if(str=="") continue;
		vector<string> tokens = string_tokenize(str, " \t");
		int N_token = tokens.size();
		ChrDir chrdir;
		chrdir.chr = tokens[0];
		chrdir.dir = tokens[1][0];
		ChrList.push_back(chrdir);
	}
	IN.close();
	if (param.Info) {
		cerr << "Finished reading chr list file." << endl;
	}
	return 1;
}

int WriteChrList ( vector<ChrDir> & ChrList ) {
	vector<ChrDir>::iterator citer;
	for ( citer = ChrList.begin(); citer != ChrList.end(); citer++ ) {
		cout << citer->chr << "\t" << citer->dir << endl;
	}
	return 1;
}


//map<string,string> GENOME;
//int ReadFasta ( string filename, vector<Seq> & SeqVec )
int ReadFasta ( string filename, map<string,string> & GENOME)
{
	// Infile
	istream *p_infile = &cin;
	ifstream infile;
	if (filename != "") {
		infile.open(filename.c_str());
		if(!infile) {
			cerr << "cannot open input file : " << filename << endl;
			exit(-1);
		}
		p_infile = &infile;
	}
	string chr = "";
	string seq = "";
	string str;
	while( getline(*p_infile, str) ){

		if(str[0] == '>'){
			if( chr != "" ) {
				GENOME[chr]=seq;
			}
			vector<string> tokens = string_tokenize( str.substr(1) );
			chr = tokens[0];
			seq = "";
		}else{
			seq += str;
		}
	}
	GENOME[chr]=seq;
	infile.close();
	if (param.Info) {
		cerr << "Finished reading input Fasta file." << endl;
	}
	return 1;
}

string Antisense ( string seq ) {
        string Antiseq = seq;
        ToUpperString ( seq );
        int len = seq.length();
        for ( int i = 0; i < len; i++ ) {
                switch (seq[i]) {
                case 'A': Antiseq[len-i-1] = 'T'; break;
                case 'C': Antiseq[len-i-1] = 'G'; break;
                case 'G': Antiseq[len-i-1] = 'C'; break;
                case 'T': Antiseq[len-i-1] = 'A'; break;
                case 'N': Antiseq[len-i-1] = 'N'; break;
                default : Antiseq[len-i-1] = 'n'; break;
                }
        }
        return Antiseq;
}


int FastaOrderByList( map<string,string> & GENOME, vector<ChrDir> & ChrList, string outfile) {
	// Outfile
	ofstream fasta_of(param.outfile.c_str());
	if( !fasta_of ) {
		cerr << "Open output file error:" << param.outfile << "\n";
		exit(-1);
	}
	vector<ChrDir>::iterator citer;
	for ( citer = ChrList.begin(); citer != ChrList.end(); citer++ ) {
		string chr = citer->chr;
		char dir = citer->dir;
		map<string,string>::iterator miter = GENOME.find(chr);
		if (miter != GENOME.end() ) {
			fasta_of << ">" << chr << endl;
			if ( citer->dir == '+' ) {
				OutputSequence( fasta_of, GENOME[chr], param.len );
			} else {
				string seq = Antisense(GENOME[chr]);
				OutputSequence( fasta_of, seq, param.len );
			}
			if ( param.Info ) {
				cerr << "Record chr: " << chr << "\t" << dir << endl;
			}
		} else {
			cerr << "[warning] The chr_id \"" << chr << "\" is not found." << endl;
		}
	}
	fasta_of.close();
	return 1;
}


void parse_command_line(int argc, char **argv)
{
	int i;
	for(i=2;i<=argc;i++)
	{
		if(argv[i-1][0] != '-') break;
		switch(argv[i-1][1])
		{
			case 'o':
				param.outfile = string(argv[i]);
				if(++i>argc) exit_with_help();
				break;
			case 'i':
				param.infile = string(argv[i]);
				if(++i>argc) exit_with_help();
				break;
			case 'c':
				param.chrlist = string(argv[i]);
				if(++i>argc) exit_with_help();
				break;
			case 'l':
				param.len= atoi(argv[i]);
				if(++i>argc) exit_with_help();
				break;
			case 'I':
				param.Info = true;
				break;
			case 'h':
				exit_with_help();
				break;
			default:
				fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
		}
	}
	if ( !param.outfile.length() ){
		exit_with_help();
	}
}



void init(){
	// default values
	param.infile = "";
	param.outfile = "";
	param.chrlist = "";
	param.len	= 50;
	param.Info = false;
}

int main(int argc, char* argv[])
{
	init();
#ifdef _DEBUG
	param.infile = string("input.fa");
	param.outfile = string("output.txt");
	param.chrlist	= string("chrlist.txt");;
	param.len	= 60;
#else
	parse_command_line(argc,argv);
#endif
	//GetGenome( param.infile.c_str() );

	map<string,string> GENOME;
	ReadFasta( param.infile, GENOME );

	vector<ChrDir> ChrList;
	ReadChrList( param.chrlist, ChrList );

	FastaOrderByList( GENOME, ChrList, param.outfile);

	if (param.Info) {
		cerr << "DONE." << endl;
	}

	return 1;
}

