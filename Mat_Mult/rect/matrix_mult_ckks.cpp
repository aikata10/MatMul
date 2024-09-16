#include "matrix_mult_ckks.h"

const int l = 8;
const int d = 8;
const int t = 4;
const int d_log =(int)(log2(d));
const int n=d*d;
int RingDim=1<<16;


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

    int rot_keys_counter= 2*d_log+2*(d)+3*(int)(log2(d));
    std::vector<int> rot_keys(rot_keys_counter);
    int ctr=0;
    for(int i=0;i<d_log;i++){
        rot_keys[ctr++]=-pow(2,i);    
        rot_keys[ctr++]=-pow(2,i)*d;     
    }
    for(int i=0;i<d;i++){
        rot_keys[ctr++]=(d)*(i+1);
        rot_keys[ctr++]=(d)*d*(i+1);
    }
    for(int i=0;i<(int)(log2(d));i++){
        rot_keys[ctr++]=-d*d*pow(2,i)+pow(2,i);
        rot_keys[ctr++]=-d*d*pow(2,i)+(pow(2,i)*t);
        rot_keys[ctr++]=n*pow(2,i);
    }

    
    m_cc->EvalRotateKeyGen(keyPair.secretKey,rot_keys);

    std::vector<double> input_matA(l*d);
    std::vector<double> input_matB(d*t);
    
    double a[l][d], b[d][t], mult[l][t];


    for(int i=0;i<(l*d);i++){
        a[(int)(i/d)][i%d]= static_cast <double> (rand()) / ( static_cast <double> (RAND_MAX/(4.4-2.2)))-1;
    	input_matA[i]=a[(int)(i/d)][i%d];
    }
    for(int i=0;i<(d*t);i++){
        b[(int)(i/t)][i%t]= static_cast <double> (rand()) / ( static_cast <double> (RAND_MAX/(4.4-2.2)))-1;
    	input_matB[i]=b[(int)(i/t)][i%t];
    }

    // Initializing elements of matrix mult to 0.
    for(int i = 0; i < l; ++i)
        for(int j = 0; j < t; ++j)
        {
            mult[i][j]=0;
        }
    

    // Multiplying matrix a and b and storing in array mult.
    //std::cout << std::endl << "Multiplying Matrix: " << std::endl;
    for(int i = 0; i < l; ++i)
        for(int j = 0; j < t; ++j)
            for(int k = 0; k < d; ++k)
            {
                mult[i][j] += a[i][k] * b[k][j];
            }


    

    size_t encodedLength = n;
    Plaintext temp;
    std::vector<std::complex<double>> finalResult2;
    
    //Matrix encoding and encryption
    
    Ciphertext<DCRTPoly> temp1, temp2;
    Plaintext plaintext  = m_cc->MakeCKKSPackedPlaintext(input_matA);
    temp1           = m_cc->Encrypt(keyPair.publicKey, plaintext);
    
    plaintext       = m_cc->MakeCKKSPackedPlaintext(input_matB);
    temp2      = m_cc->Encrypt(keyPair.publicKey, plaintext);
    
    
    
    //---------------------- Multiplication evaluation -----------------------------------
    
   // std::cout << std::endl << "EVAL MULT STARTS " << std::endl;
   clock_t begin = clock();


    std::vector<double> mask1(d*d*d,0);
    std::vector<std::vector<double>> mask;
    mask.push_back(mask1);
    mask.push_back(mask1);
    
    
    Ciphertext<DCRTPoly> out;
    std::vector<Plaintext> plaintext_mask(2);
    
    for(int kloop=0;kloop<2*d;kloop++){   

            int tem= l*(1-(kloop%2))+t*(kloop%2);     
		        for(int i=0;i<tem;i++){
		                int t=(kloop%2)*(i)+(1-(kloop%2))*(i*d);
			    	mask[(kloop%2)][t+(int)(kloop/2)*d*d]=1 ;
			    }
		    }

    
    
    #pragma omp parallel for
    for(int j=0;j<2;j++){
    	plaintext_mask[j]  = m_cc->MakeCKKSPackedPlaintext(mask[j]);
    }
    
    for(int n_runs=0;n_runs<100;n_runs++){
    	m_MatrixAC=temp1;
    m_MatrixBC=temp2;
    for(int k=0;k<d_log;k++)
        {   
             m_MatrixAC=m_cc->EvalAdd(m_MatrixAC,m_cc->EvalRotate(m_MatrixAC,-d*d*pow(2,k)+pow(2,k)));
            m_MatrixBC=m_cc->EvalAdd(m_MatrixBC,m_cc->EvalRotate(m_MatrixBC,-d*d*pow(2,k)+(pow(2,k)*t)));
        }

    
    

    //std::cout << std::endl << " Initial A,B packing done" << std::endl;
    	
    
    
	
    std::vector<Ciphertext<DCRTPoly>> ab1(2),ab2(2);
    ab1[0]=m_MatrixAC;
    ab1[1]=m_MatrixBC;

    //std::cout << std::endl << " LOOP initial Rotates done " << std::endl;
	
             
    //#pragma omp parallel for
	 int j=0;
	    	ab1[j]  =m_cc->EvalMult(ab1[j],plaintext_mask[j]);
	    	m_cc->ModReduceInPlace(ab1[j]);
            //std::cout << std::endl << " LOOP initial Rotates done "<< j << " " << std::endl;
	
	    	for(int kl=0;kl<int(log2(d));kl++){
		    	int mo=-1*pow(2,kl);
		    	ab2[j]=m_cc->EvalRotate(ab1[j],mo*pow(d,j));
			    ab1[j]=m_cc->EvalAdd(ab1[j],ab2[j]);	 
		    }

    j=1;
    ab1[j]  =m_cc->EvalMult(ab1[j],plaintext_mask[j]);
	    	m_cc->ModReduceInPlace(ab1[j]);
           // std::cout << std::endl << " LOOP initial Rotates done "<< j << " " << std::endl;
	
	    	for(int kl=0;kl<int(log2(d));kl++){
		    	int mo=-1*pow(2,kl);
		    	ab2[j]=m_cc->EvalRotate(ab1[j],mo*pow(d,j));
			    ab1[j]=m_cc->EvalAdd(ab1[j],ab2[j]);	 
		    }

	
	    
    out=m_cc->EvalMult(ab1[0],ab1[1]);
    m_cc->ModReduceInPlace(out);

    for(int k=0;k<(int)(log2(d));k++)
        {   
             out=m_cc->EvalAdd(out,m_cc->EvalRotate(out,pow(2,k)*n)); 
        } 
        
       }

    
    m_OutputC=out;
    
    clock_t end = clock();
   double time_spent = (double)(end - begin) / (100*CLOCKS_PER_SEC);
    


    
    
    
    
    //---------------------- Multiplication evaluation done -----------------------------------
       
    
    //Result decryption
    
    Plaintext plaintextDec;
    m_cc->Decrypt(keyPair.secretKey, m_OutputC, &plaintextDec);
    plaintextDec->SetLength(encodedLength);



    std::vector<std::complex<double>> finalResult = plaintextDec->GetCKKSPackedValue();
    std::complex<double> error=0;
    for(int i=0;i<(l*t);i++){
    	error+=(finalResult[i+((int)(i/t)*(int)(d-t))])-(std::complex<double>)(mult[(int)(i/t)][i%t]);
    	
    	}
    
    //std::cout << "Actual output\n\t" << finalResult << std::endl << std::endl;
    
     std::cout << "Error\n\t" << error << std::endl << std::endl;
     std::cout << "Time\n\t" << time_spent << std::endl << std::endl;
    
    
    
}

void MatrixMultCKKS::deserializeOutput()
{
    

}
