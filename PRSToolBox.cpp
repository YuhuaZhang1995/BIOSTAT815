#include <iostream>
#include <stdlib.h>
#include <ctpl_stl.h>
#include <sstream>
#include <string>
#include <fstream>
#include <RcppArmadillo.h>
#include <queue>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

int get_snpnum(const char* filename){
	//This function is to get the number of SNPs

	/*Read vcf file line by line*/
	std::ifstream fin(filename,std::ios_base::in);
	std::string line;
	int snpnum=0;
	while (std::getline(fin,line)){
		/*Remove the header of the .vcf file*/
		if (line.find("##")==std::string::npos){
			
			/*get the number of SNPs*/
			if(line.find("#")==std::string::npos){
				snpnum++;
			}
		}
	}
	fin.close();
	return(snpnum);
}

int get_ind(const char* filename){
	//This function is to get the number of individuals

	/*Read vcf file line by line*/
	std::ifstream fin(filename,std::ios_base::in);
	std::string line;
	int count=0;
	while (std::getline(fin,line)){
		/*Remove the header of the .vcf file*/
		std::istringstream iss(line);
		std::string word;
		if (line.find("##")==std::string::npos){
			
			/*First line of the vcf file, get the number of individuals*/
			if (line.find("#")!=std::string::npos){
				while(iss>>word){
					if ((word!="#CHROM")&(word!="POS")&(word!="ID")&(word!="REF")&(word!="ALT")&(word!="QUAL")&(word!="FILTER")&(word!="INFO")&(word!="FORMAT")){
						count++;
					}
				}
			}
			break;
		}
	}
	fin.close();
	return(count);
}

arma::mat parse_file(const char* filename, int start_line, int end_line, std::queue<String>& id_list){
	//This function is to read the file and return the parsed genotype matrix
	/*Get the number of individuals and number of SNPs*/
	//int snpnum=get_snpnum(filename);
	int indnum=get_ind(filename);

	std::ifstream fin(filename,std::ios_base::in);
	std::string line;
	arma::mat geno((end_line-start_line+1),indnum);

	int count=0;
	while (std::getline(fin,line)){
		std::istringstream iss(line);
		std::string word;
		/*Remove the header of the .vcf file*/
		if (line.find("##")==std::string::npos){
			if(line.find("#")==std::string::npos){
				if((start_line<=count)&(count<(end_line+1))){
				/*store the parsed matrix into geno mat*/
				int tmp_count=0;
				int tmp_count2=0;
				while(iss>>word){
					if (word.find("/")!=std::string::npos){
						std::size_t pos=word.find("/");
						std::string tmp_str=word.substr(0,pos);
						try{geno((count-start_line),tmp_count)=stod(tmp_str);}
						catch(std::exception& e){geno((count-start_line),tmp_count)=0.0;}
						tmp_count++;
					}
					if (tmp_count2==2){
						id_list.push(word);
					}
					tmp_count2+=1;
				}
				}count++; //Rcout<<count<<std::endl;
			}
		}
	}
	fin.close();
	return(geno);
}

int get_weighted_snpnum(const char* filename){
	//This function is to get the number of traits in weighted file

	int count=0;
	std::ifstream fin2(filename,std::ios_base::in);
	std::string line;
	while (std::getline(fin2,line)){
		std::istringstream iss(line);
		std::string word;
		//Rcout<<line<<std::endl;
		while(iss>>word){
			count++;	
		}
		break;
	}
	int tmp=count-1;
	fin2.close();
	return(tmp);
}

arma::mat weights_file(const char* filename, std::queue<String>& id_list){
	//This function is to read weights file and return it as the arma::mat
	//The returned weight matrix only contains those variants in the id_list

	int num_weights=get_weighted_snpnum(filename); //number of traits
	int snpnum=id_list.size(); //number of snps

	/*Generate the matrix*/
	arma::mat weights(snpnum,num_weights);
	std::string line;
	std::string str_in_id_list;
	//Rcout<<"label"<<std::endl;
	int count2=0; //row index in weight matrix
	while (id_list.size()>0){
		str_in_id_list=id_list.front();
		id_list.pop();
		int count=0;
		std::ifstream fin(filename,std::ios_base::in);
		while (std::getline(fin,line)){
			std::istringstream iss(line);
			std::string word;
			//std::string str_in_id_list;
			if(count>0){
				int tmp_count=0;
				bool FLAG=false;
				while(iss>>word){
					if (tmp_count==0){
						if (word==str_in_id_list){
							//Rcout<<word<<str_in_id_list<<std::endl;
							FLAG=true;
							count2++;
							//Rcout<<"label"<<std::endl;
						}
					}
					if(FLAG==false){break;}
					if ((tmp_count>0) & (FLAG==true)){
						//Rcout<<"label"<<std::endl;
						weights(count2-1,tmp_count-1)=stod(word);
						//Rcout<<count2-1<<std::endl;
					}tmp_count++;
				}
			}
			count++;
		}fin.close();
	}
	return(weights);
}

//' An Rcpp implementation of calculating polygenic risk score.
//'
//' @param geno_file the genotype file (in .vcf format) that contains individual level genotype information
//' @param weight_file the weight file, which is a variant by individual matrix
//' @param thread_pool_num, number of thread pools
//' @param num_thread number of chunks to separate the genotype files into
//' @return An individual by trait matrix containing the PRS.
// [[Rcpp::export]]
arma::mat cal_PRS(const char* geno_file, const char* weight_file,int thread_pool_num, int num_thread){
	/*Get the number of SNPs*/
	int snpnum=get_snpnum(geno_file);
	int indnum=get_ind(geno_file);
	//arma::mat weights(snpnum,1);
	//weights=weights_file(weight_file);

	/*parse the files based on the number of SNPS*/
	int mat_size=snpnum/num_thread; //number of the SNPs to be enrolled in on matrix
	int trait=get_weighted_snpnum(weight_file); //number of traits
	ctpl::thread_pool p(thread_pool_num); //number of thread pools
	arma::mat tmp_result(indnum,trait); //store the final result
	tmp_result.zeros();
	for (int i=0;i<num_thread;i++){
		if(i<num_thread-1){
	/*Get the genotype matrix*/
		std::queue<String> id_list;
		arma::mat geno=parse_file(geno_file,(i*mat_size),(i*mat_size+mat_size-1),id_list);
		//std::string tmp_word;
		//tmp_word=id_list.front();
		//Rcout<<tmp_word<<std::endl;		
	/*Get the weights*/
		arma::mat weight(id_list.size(),trait);
		weight=weights_file(weight_file,id_list);
	/*Calculate PRS*/
		p.push([&tmp_result,geno,weight](int){
			tmp_result+=(geno.t()*weight);
		});}
	/*Deal with the last few lines of the genotype files*/
		if(i==num_thread-1){
			std::queue<String> id_list;
			arma::mat geno=parse_file(geno_file,(i*mat_size),(snpnum-1),id_list);
		
			arma::mat weight(id_list.size(),trait);
			weight=weights_file(weight_file,id_list);
			/*Calculate PRS*/
			p.push([&tmp_result,geno,weight](int){
				tmp_result+=(geno.t()*weight);
			});
		}
	}
	//Rcout<<tmp_result<<std::endl;
	return tmp_result;
}
