#include "matrix_mult_ckks.h"


const int d = 128;
const int d_sq = (int)(sqrt(d));
const int d_log =(int)(log2(d));
const int n=d*d;
int RingDim=1<<16;
int packing=std::min(d,(int)(RingDim/(2*d*d)));
int n_ct=(int)(d/packing);


MatrixMultCKKS::MatrixMultCKKS(std::string ccLocation, std::string pubKeyLocation, std::string multKeyLocation,
                               std::string rotKeyLocation,
                               std::string matrixALocation,
                               std::string matrixBLocation,
                               std::string outputLocation) : m_PubKeyLocation(pubKeyLocation),
                                                             m_MultKeyLocation(multKeyLocation),
                                                             m_RotKeyLocation(rotKeyLocation),
                                                             m_CCLocation(ccLocation),
                                                             m_MatrixALocation(matrixALocation),
                                                             m_MatrixBLocation(matrixBLocation),
                                                             m_OutputLocation(outputLocation)
{

    //initCC();
};

void MatrixMultCKKS::initCC()
{
    
}

void MatrixMultCKKS::eval()
{
    // You should implement this function
    
    
    
    //generate keys here and serialize them into file.... the  uncomment the init in the top function to load keys from the file always.
    std::cout << "--------------------------------- EVAL Matrix Multiplication ---------------------------------"
              << std::endl;
    CCParams<CryptoContextCKKSRNS> parameters;

    // We set a smaller ring dimension to improve performance for this example.
    // In production environments, the security level should be set to
    // HEStd_128_classic, HEStd_192_classic, or HEStd_256_classic for 128-bit, 192-bit,
    // or 256-bit security, respectively.
    parameters.SetSecurityLevel(HEStd_128_classic);
    
    parameters.SetRingDim(RingDim);
    //parameters.SetBatchSize(n);
    usint scalingModSize = 50;
    //usint firstModSize   = 60;

    parameters.SetScalingModSize(scalingModSize);
    //parameters.SetFirstModSize(firstModSize);

    // The multiplicative depth depends on the polynomial degree.
    uint32_t multDepth = 2;

    parameters.SetMultiplicativeDepth(multDepth);
    CryptoContext<DCRTPoly> m_cc = GenCryptoContext(parameters);
    m_cc->Enable(PKE);
    m_cc->Enable(KEYSWITCH);
    m_cc->Enable(LEVELEDSHE);

    KeyPair<DCRTPoly> keyPair = m_cc->KeyGen();
    // We need to generate mult keys
    m_cc->EvalMultKeyGen(keyPair.secretKey);

    int rot_keys_counter= 2*d_log+2*(n_ct-1)+3*(int)(log2(packing));

    std::cout << std::endl << "Total ROTATIONS:" << ( n_ct*2*d_log+2*(n_ct-1)+3*(int)(log2(packing))) << std::endl;
    std::cout << std::endl << "Total MULTIPLICATIONS:" << ( 3*n_ct) << std::endl;
    std::cout << std::endl << "Total PLAINTEXT MULTIPLICATIONS:" << ( 2*n_ct) << std::endl;
    std::cout << std::endl << "Total CIPHERTEXT MULTIPLICATIONS:" << ( n_ct) << std::endl;

    std::vector<int> rot_keys(rot_keys_counter);

    for(int i=0;i<d_log;i++){
        rot_keys[2*i]=-pow(2,i);
        rot_keys[2*i+1]=-pow(2,i)*d;
    }
    for(int i=0;i<n_ct-1;i++){
        rot_keys[2*i+2*d_log]=(packing)*(i+1);
        rot_keys[2*i+2*d_log+1]=(packing)*d*(i+1);
    }
    for(int i=0;i<(int)(log2(packing));i++){
        rot_keys[3*i+2*d_log+2*n_ct-2]=-n*pow(2,i)+pow(2,i);
        rot_keys[3*i+2*d_log+2*n_ct-2+1]=-n*pow(2,i)+pow(2,i)*d;
        rot_keys[3*i+2*d_log+2*n_ct-2+2]=n*pow(2,i);
    }

    
    


    //std::cout  << "rotkey gen: "  << rot_keys << std::endl;
    m_cc->EvalRotateKeyGen(keyPair.secretKey,rot_keys);

    std::vector<double> input_matA(n);
    std::vector<double> input_matB(n);
    
    double a[d][d], b[d][d], mult[d][d];

    //std::cout << std::endl << " Matrix gen: " << std::endl;
    
    for(int i=0;i<(n);i++){
    	a[(int)(i/d)][i%d]= static_cast <double> (rand()) / ( static_cast <double> (RAND_MAX/(4.4-2.2)))-1;
    	input_matA[i]=a[(int)(i/d)][i%d];
    	b[(int)(i/d)][i%d]= static_cast <double> (rand()) / ( static_cast <double> (RAND_MAX/(4.4-2.2)))-1;
    	input_matB[i]=b[(int)(i/d)][i%d];
    }
    //std::cout << "Matrix 1 \n\t" << input_matA << std::endl << std::endl;
    //std::cout << "Matrix 2 \n\t" << input_matB << std::endl << std::endl;
    
    // Initializing elements of matrix mult to 0.
    for(int i = 0; i < d; ++i)
        for(int j = 0; j < d; ++j)
        {
            mult[i][j]=0;
        }
    

    // Multiplying matrix a and b and storing in array mult.
    //std::cout << std::endl << "Multiplying Matrix: " << std::endl;
    for(int i = 0; i < d; ++i)
        for(int j = 0; j < d; ++j)
            for(int k = 0; k < d; ++k)
            {
                mult[i][j] += a[i][k] * b[k][j];
            }

    // Displaying the multiplication of two matrix.
    
    //std::cout << std::endl << "Output Matrix: " << std::endl;
    /*for(int i = 0; i < d; ++i)
    for(int j = 0; j < d; ++j)
    {
        std::cout << " " << mult[i][j];
        if(j == (d-1))
        
            std::cout << std::endl;
    }*/
    
    
    

    size_t encodedLength = input_matA.size();
    
    //Matrix encoding and encryption
    Ciphertext<DCRTPoly> temp1, temp2;
    Plaintext plaintext  = m_cc->MakeCKKSPackedPlaintext(input_matA);
    temp1           = m_cc->Encrypt(keyPair.publicKey, plaintext);
    
    plaintext       = m_cc->MakeCKKSPackedPlaintext(input_matB);
    temp2      = m_cc->Encrypt(keyPair.publicKey, plaintext);
    
    
    
    //---------------------- Multiplication evaluation -----------------------------------
    
    std::cout << std::endl << "EVAL MULT STARTS " << std::endl;
    
    clock_t begin = clock();


    std::vector<double> mask1(packing*n,0);
    std::vector<std::vector<double>> mask;
    mask.push_back(mask1);
    mask.push_back(mask1);
    
    
    std::vector<Ciphertext<DCRTPoly>> out(n_ct);
    std::vector<Plaintext> plaintext_mask(2);
    #pragma omp parallel for collapse(2)
    for(int k=0;k<2*packing;k++){
		        for(int i=0;i<d;i++){
		                int t=(k%2)*(i)+(1-(k%2))*(i*d);
			    	mask[(k%2)][t+(int)(k/2)*n]=1 ;
			    }
		    }
    
    
    #pragma omp parallel for
    for(int j=0;j<2;j++){
    	plaintext_mask[j]  = m_cc->MakeCKKSPackedPlaintext(mask[j]);
    }
    
    for(int nruns=0;nruns<100;nruns++){
    
    m_MatrixAC=temp1;
    m_MatrixBC=temp2;
   
    
   
    for(int k=0;k<(int)(log2(packing));k++)
    {   
	    m_MatrixAC=m_cc->EvalAdd(m_MatrixAC,m_cc->EvalRotate(m_MatrixAC,-n*pow(2,k)+pow(2,k)));
	    m_MatrixBC=m_cc->EvalAdd(m_MatrixBC,m_cc->EvalRotate(m_MatrixBC,-n*pow(2,k)+pow(2,k)*d));
    }
    
    

    //std::cout << std::endl << " Initial A,B packing done" << std::endl;
    	
    
    
	
     #pragma omp parallel for shared(m_MatrixAC,m_MatrixBC,plaintext_mask)
    for(int t=0;t<n_ct;t++){
        std::vector<Ciphertext<DCRTPoly>> ab1(2),ab2(2);
	    if(t!=0){    
		        #pragma omp parallel sections
			    {   
				    #pragma omp section
				   ab1[0]=m_cc->EvalRotate(m_MatrixAC,packing*t);
				    #pragma omp section
				   ab1[1]=m_cc->EvalRotate(m_MatrixBC,packing*t*d);
			    }        
	     }
	    else{
		    #pragma omp parallel sections
			    {   
				    #pragma omp section
				   ab1[0]=m_MatrixAC;
				    #pragma omp section
				   ab1[1]=m_MatrixBC;
			    }	        
	     }

        // std::cout << std::endl << " LOOP initial Rotates done " << std::endl;
	
             
          //#pragma omp parallel for
	 for(int j=0;j<2;j++){
	    	ab1[j]  =m_cc->EvalMult(ab1[j],plaintext_mask[j]);
	    	m_cc->ModReduceInPlace(ab1[j]);
	    	for(int k=0;k<d_log;k++){
		    	int l=-1*pow(2,k);
		    	ab2[j]=m_cc->EvalRotate(ab1[j],l*pow(d,j));
			ab1[j]=m_cc->EvalAdd(ab1[j],ab2[j]);	 
		    }
	 }
	      
	    
	     out[t]=m_cc->EvalMult(ab1[0],ab1[1]);
    	    m_cc->ModReduceInPlace(out[t]);
	    
	   
    }
    
			   
	    for(int i=1;i<=(int)(log2(n_ct));i++){
            #pragma omp parallel for 
	        for(int t=0;t<(int)(n_ct/pow(2,i));t++){
                m_cc->EvalAddInPlace(out[t], out[t+(int)(n_ct/pow(2,i))]);
            }
        }


    for(int k=0;k<(int)(log2(packing));k++)
        {   
             out[0]=m_cc->EvalAdd(out[0],m_cc->EvalRotate(out[0],pow(2,k)*n)); 
        }  

    m_OutputC=out[0];
    
    }
    
    clock_t end = clock();
    double time_spent = (double)(end - begin) / (100*CLOCKS_PER_SEC);
    
    
    
    

     
   
	   /* Plaintext temp;
	   m_cc->Decrypt(keyPair.secretKey, m_OutputC, &temp);
	    temp->SetLength(encodedLength);
    	    std::vector<std::complex<double>> finalResult2 = temp->GetCKKSPackedValue();   
           std::cout << " Input 1 \n\t" << finalResult2 << std::endl << std::endl;
           //m_cc->Decrypt(keyPair.secretKey, b1, &temp);
	   */
	    
	    
	    
	   
	    
    

    

   // std::cout << " NoiseScale output \n\t" << m_OutputC->GetNoiseScaleDeg() << std::endl << std::endl;
    



    
    
    
    
    //---------------------- Multiplication evaluation done -----------------------------------
       
    
    //Result decryption
    
    Plaintext plaintextDec;
    m_cc->Decrypt(keyPair.secretKey, m_OutputC, &plaintextDec);
    plaintextDec->SetLength(encodedLength);



    std::vector<std::complex<double>> finalResult = plaintextDec->GetCKKSPackedValue();
    std::complex<double> error=0;
    for(int i=0;i<(n);i++){
    	error+=(finalResult[i])-(std::complex<double>)(mult[(int)(i/d)][i%d]);
    	
    	}
    
    //std::cout << "Actual output\n\t" << finalResult << std::endl << std::endl;
    
     std::cout << "Error\n\t" << error << std::endl << std::endl;
     std::cout << "Time\n\t" << time_spent << std::endl << std::endl;
    
    
    
}

void MatrixMultCKKS::deserializeOutput()
{

}
