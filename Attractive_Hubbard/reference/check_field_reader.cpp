#include <iostream>
#include <fstream>
int main(int argc,char**argv){std::ifstream in(argv[1]);double x;int n=0;while(in>>x)++n;std::cout<<"read_count="<<n<<" fail="<<in.fail()<<" eof="<<in.eof()<<"\n";}
