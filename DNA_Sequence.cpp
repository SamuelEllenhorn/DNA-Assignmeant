#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
using namespace std;
/*
Samuel Ellenhorn
2295046
Ellenhorn@chapman.edu
Computer Science 350
Programming Assignment 1
*/

//Global Variables declared here so they can bet set later at their respective times.
int numlines = 0;
long sumSquares = 0;
long numChar = 0;
double variance = 0;
double mean = 0;
double standardDeviation = 0;

double probA = 0.0;
double probT = 0.0;
double probC = 0.0;
double probG = 0.0;

int Acount = 0;
int Tcount = 0;
int Ccount = 0;
int Gcount = 0;

int AA= 0;
int AC= 0;
int AT= 0;
int AG = 0;
int TA= 0;
int TC= 0;
int TT= 0;
int TG = 0;
int CA= 0;
int CC= 0;
int CT= 0;
int CG = 0;
int GA= 0;
int GC= 0;
int GT= 0;
int GG = 0;

double probAA = 0.0;
double probAC = 0.0;
double probAT = 0.0;
double probAG = 0.0;
double probTA = 0.0;
double probTC = 0.0;
double probTT = 0.0;
double probTG = 0.0;
double probCA = 0.0;
double probCC = 0.0;
double probCT = 0.0;
double probCG = 0.0;
double probGA = 0.0;
double probGC = 0.0;
double probGT = 0.0;
double probGG = 0.0;


//The following boolean defenitions funcition is used to void the capitilization state from userInput.
bool is_A(char c) {
    return c == 'a' || c == 'A';
}
bool is_T(char c) {
    return c == 't' || c == 'T';
}
bool is_C(char c) {
    return c == 'c' || c == 'C';
}
bool is_G(char c) {
    return c == 'g' || c == 'G';
}


//The purpose of this function is to set all global var to 0.
void reset(){
    numlines = 0;
    sumSquares = 0;
    numChar = 0;
    variance = 0;
    mean = 0;
    standardDeviation = 0.0;

    Acount = 0;
    Tcount = 0;
    Ccount = 0;
    Gcount = 0;

    AA = 0;
    AC = 0;
    AT = 0;
    AG = 0;
    TA = 0;
    TC = 0;
    TT = 0;
    TG = 0;
    CA = 0;
    CC = 0;
    CT = 0;
    CG = 0;
    GA = 0;
    GC = 0;
    GT = 0;
    GG = 0;

    probAA = 0.0;
    probAC = 0.0;
    probAT = 0.0;
    probAG = 0.0;
    probTA = 0.0;
    probTC = 0.0;
    probTT = 0.0;
    probTG = 0.0;
    probCA = 0.0;
    probCC = 0.0;
    probCT = 0.0;
    probCG = 0.0;
    probGA = 0.0;
    probGC = 0.0;
    probGT = 0.0;
    probGG = 0.0;

}


//The purpose of this function is to analyze the file and to set global variables.
void process_line(string s) {
    numChar += s.length();
    sumSquares += s.length() * s.length();
    for (size_t i = 0; i < s.length(); i++) {
        if (is_A(s[i]))
            Acount++;
        if (is_T(s[i]))
            Tcount++;
        if (is_C(s[i]))
            Ccount++;
        if (is_G(s[i]))
            Gcount++;
    //The following code is used to calculate the nucleotide bigram variables
        if (s.length() > 1) {
            for(size_t i = 1; i < s.length(); i++) {
                if(is_A(s[i-1]) && is_A(s[i]))
                    AA++;
                if(is_A(s[i-1]) && is_C(s[i]))
                    AC++;
                if(is_A(s[i-1]) && is_T(s[i]))
                    AT++;
                if(is_A(s[i-1]) && is_G(s[i]))
                    AG++;
    ////////////////////////////////////////
                if(is_T(s[i-1]) && is_A(s[i]))
                    TA++;
                if(is_T(s[i-1]) && is_C(s[i]))
                    TC++;
                if(is_T(s[i-1]) && is_T(s[i]))
                    TT++;
                if(is_T(s[i-1]) && is_G(s[i]))
                    TG++;
    /////////////////////////////////////////
                if(is_C(s[i-1]) && is_A(s[i]))
                    CA++;
                if(is_C(s[i-1]) && is_C(s[i]))
                    CC++;
                if(is_C(s[i-1]) && is_T(s[i]))
                    CT++;
                if(is_C(s[i-1]) && is_G(s[i]))
                    CG++;
    ////////////////////////////////////////
                if(is_G(s[i-1]) && is_A(s[i]))
                    GA++;
                if(is_G(s[i-1]) && is_C(s[i]))
                    GC++;
                if(is_G(s[i-1]) && is_T(s[i]))
                    GT++;
                if(is_G(s[i-1]) && is_G(s[i]))
                    GG++;
        }
    }
    numlines++;
}


//The purpose of this function is to compute the variable values which require the entire data set to be processed prior.
void compute_stats() {
    mean = (double)numChar/(double)numlines;




    variance = static_cast<double>(sumSquares) / (double)numlines - mean * mean;
    // This is part Two of the standardDeviation calculation.
    standardDeviation = sqrt(variance);


    //computing nucleotide probability
    probA = Acount/(double)numChar;
    probT = Tcount/(double)numChar;
    probC = Ccount/(double)numChar;
    probG = Gcount/(double)numChar;

    //computing nucleotide-bigram probability
    double numBiGram = numChar - numlines;
    probAA = AA/numBiGram;
    probAC = AC/numBiGram;
    probAT = AT/numBiGram;
    probAG = AG/numBiGram;
    probTA = TA/numBiGram;
    probTC = TC/numBiGram;
    probTT = TT/numBiGram;
    probTG = TG/numBiGram;
    probCA = CA/numBiGram;
    probCC = CC/numBiGram;
    probCT = CT/numBiGram;
    probCG = CG/numBiGram;
    probGA = GA/numBiGram;
    probGC = GC/numBiGram;
    probGT = GT/numBiGram;
    probGG = GG/numBiGram;

}

//C = sqrt(-2 ln (a)) * cos(2πb)

//D = σC + μ

char random_nucleotide() {
    double r = rand()*1.0/RAND_MAX;

}

//the following code is used to generate the new DNA strings while factoring in previoes nucleotide/bigram probability.
char random_nucleotide(char prev) {
    double r = rand()*1.0/RAND_MAX;
    if(prev == 'A') {
        double s = probAA + probAT + probAC + probAG;
        s = r*s;
        if (s <= probAA)
            return 'A';
        else if (s <= probAT + probAA)
            return 'T';
        else if (s <=  probAA + probAT + probAC)
            return 'C';
        else
            return 'G';
    }
    else if(prev == 'T') {
        double s = probTA + probTT + probTC + probTG;
        s = r*s;
        if (s <= probTA)
            return 'A';
        else if (s <= probTT + probTA)
            return 'T';
        else if (s <=  probTA + probTT + probTC)
            return 'C';
        else
            return 'G';
    }
    else if(prev == 'C') {
        double s = probCA + probCT + probCC + probCG;
        s = r*s;
        if (s <= probCA)
            return 'A';
        else if (s <= probCT + probCA)
            return 'T';
        else if (s <=  probCA + probCT + probCC)
            return 'C';
        else
            return 'G';
    }
    else if(prev == 'G') {
        double s = probGA + probGT + probGC + probGG;
        s = r*s;
        if (s <= probGA)
            return 'A';
        else if (s <= probGT + probGA)
            return 'T';
        else if (s <=  probGA + probGT + probGC)
            return 'C';
        else
            return 'G';
    } else {
        if (r <= probA)
            return 'A';
        else if (r <= probT + probA)
            return 'T';
        else if (r <= probT + probC + probA)
            return 'C';
        else
            return 'G';
    }
}

string generate_string(){
    double a = rand()*1.0/RAND_MAX;
    double b = rand()*1.0/RAND_MAX;

    double c = sqrt(-2 *log (a)) * cos(2*M_PI*b);
    int d = standardDeviation * c + mean;
    // this is a random length normally distributed

    string s = "";
    char prev = ' ';
    for(int i = 0; i < d; i++) {
        prev = random_nucleotide(prev);
        s += prev;
    }

    return s;

}


//This function will print my name and student I
//& : pass by refference so that modications i make here will effect the fout object.
void output_name_info(ofstream & fout) {
    fout << "Name: Samuel Ellenhorn" << endl;
    fout << "Student ID: " << endl;
}


void output_stats(bool is_append_name_info) {
    ofstream fout;
    fout.open("SamuelEllenhorn.out",  ofstream::app);

    if (is_append_name_info)// this is the first time opening printing to file )
        output_name_info(fout);

    fout << "mean: " << mean << endl;
    fout << "variance: " << variance << endl;
    fout << "Standard Deviation: "<< standardDeviation << endl;
    fout << "sum: " << numChar << endl;
    fout << "The relative probability of each nucletide:"<< endl;
    fout << "A: " <<Acount<< endl;
    fout << "C: " <<Ccount<< endl;
    fout << "T: " <<Tcount<< endl;
    fout << "G: " <<Gcount<< endl;
    fout << "The relative probability of each nucletide bigram:"<< endl;

    fout << "AA: " << probAA << endl;
    fout << "AT: " << probAT << endl;
    fout << "AC: " << probAC << endl;
    fout << "AG: " << probAG << endl;

    fout << "TA: " << probTA << endl;
    fout << "TT: " << probTT << endl;
    fout << "TC: " << probTC << endl;
    fout << "TG: " << probTG << endl;

    fout << "CA: " << probCA << endl;
    fout << "CT: " << probCT << endl;
    fout << "CC: " << probCC << endl;
    fout << "CG: " << probGG << endl;

    fout << "GA: " << probGA << endl;
    fout << "GT: " << probGT << endl;
    fout << "GC: " << probGC << endl;
    fout << "GG: " << probGG << endl;

    for(int i = 0; i < 1000; i++)
        fout << generate_string() << endl;


}

// returns true if successfully processed file
bool proccess_file(string fileName){
    string line;
    //opens the command line param
    ifstream fin(fileName);

    if(!fin) {
        cout << "Unable to open file" << endl;
        return false;
    }

    while(getline(fin, line)) {
        // line now contains the next lineCounter
        process_line(line);
    }

    // Done reading file. compute output_stats
    compute_stats();
    fin.close();
    return true;
}









int main(int argc, char** argv)
{
    string again = "no";
    string fileName = argv[1];
    bool is_append_name_info = true;
    do {

        if(proccess_file(fileName)) {
            output_stats(is_append_name_info);
        }

        
        cout << "Would you like to process another list?(yes/no)" << endl;
        cin >> again;
        if(again == "yes") {
            cout << "What is the name of the file?" <<endl;
            cin >> fileName;
        }

        is_append_name_info = false; // only append name and id at the top
        reset();
    } while(again == "yes");


}
