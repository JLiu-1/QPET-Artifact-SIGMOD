#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "QoZ/api/sz.hpp"

#define SZ_FLOAT 0
#define SZ_DOUBLE 1
#define SZ_UINT8 2
#define SZ_INT8 3
#define SZ_UINT16 4
#define SZ_INT16 5
#define SZ_UINT32 6
#define SZ_INT32 7
#define SZ_UINT64 8
#define SZ_INT64 9

void usage() {
    printf("Note: HPEZ command line arguments are backward compatible with SZ2/3, \n");
    //printf("      use -h2 to show the supported SZ2 command line arguments. \n");
    printf("Usage: hpez <options>\n");
    printf("Options:\n");
    printf("* general options:\n");
    printf("	-h: print the help information\n");
    printf("	-h2: print the help information for SZ2 style command line\n");
    printf("	-v: print the version number\n");
    printf("	-a : print compression results such as distortions\n");
    printf("* input and output:\n");
    printf("	-i <path> : original 3 binary input files.\n");
    printf("	-o <path> : 3 compressed output files, default in binary format\n");
    printf("	-z <path> : 3 compressed output (w -i) or input (w/o -i) files\n");
    printf("	-t : store compressed output file in text format\n");
//    printf("	-p: print meta data (configuration info)\n");
    printf("* data type:\n");
    printf("	-f: single precision (float type)\n");
    printf("	-d: double precision (double type)\n");
    printf("	-I <width>: integer type (width = 32 or 64)\n");
    printf("* configuration file: \n");
    printf("	-c <configuration file> : configuration file qoz.config\n");
    printf("* error control: (the error control parameters here will overwrite the setting in sz.config)\n");
    printf("	-M <error control mode> <error bound (optional)> \n");
    printf("	error control mode as follows: \n");
    printf("		ABS (absolute error bound)\n");
    printf("		REL (value range based error bound, so a.k.a., VR_REL)\n");
    printf("		PSNR (peak signal-to-noise ratio)\n");
    printf("		NORM (norm2 error : sqrt(sum(xi-xi')^2)\n");
    printf("		ABS_AND_REL (using min{ABS, REL})\n");
    printf("		ABS_OR_REL (using max{ABS, REL})\n");
    printf("	error bound can be set directly after the error control mode, or separately with the following options:\n");
    printf("		-A <absolute error bound>: specifying absolute error bound\n");
    printf("		-R <value_range based relative error bound>: specifying relative error bound\n");
//    printf("		-P <point-wise relative error bound>: specifying point-wise relative error bound\n");
    printf("		-S <PSNR>: specifying PSNR\n");
    printf("		-N <normErr>: specifying normErr\n");
    printf("    -q <level>: level for activation of qoz/hpez features.\n");
    printf("        Level 0: SZ3.1 (No QoZ features).\n");
    printf("        Level 1: QoZ1 (Anchor-point-based level-wise interpolation tuning).\n");
    printf("        Level 2: HPEZ l2 (level 1 plus multi-dim interp, natural cubic spline, interpolation re-ordering and dynamic dimension freezing).\n");
    printf("        Level 3: HPEZ l3 (level 2 plus dynamic dimension weights).\n");
    printf("        Level 4: HPEZ l4 (level 3 plus block-wise interpolation tuning).\n");
    printf("    -T <QoZ tuning target> \n");
    printf("    tuning targets as follows: \n");
    printf("        PSNR (peak signal-to-noise ratio)\n");
    printf("        CR (compression ratio)\n");
    printf("        SSIM (structural similarity)\n");
    printf("        AC (autocorrelation)\n");
    printf("    -C <anchor stride> : stride of anchor points.\n");
    printf("    -B <sampling block size> : block size of sampled data block for auto-tuning.\n");
    printf("    -u <quantile> : eb quantile\n");
    printf("* dimensions: \n");
    printf("	-1 <nx> : dimension for 1D data such as data[nx]\n");
    printf("	-2 <nx> <ny> : dimensions for 2D data such as data[ny][nx]\n");
    printf("	-3 <nx> <ny> <nz> : dimensions for 3D data such as data[nz][ny][nx] \n");
    printf("	-4 <nx> <ny> <nz> <np>: dimensions for 4D data such as data[np][nz][ny][nx] \n");
    printf("* examples: \n");
    printf("	hpez -f -i test.dat    -z test.dat.qoz     -3 8 8 128 -M ABS 1e-3 -q 3\n");
    printf("	hpez -f -z test.dat.qoz -o test.dat.qoz.out -3 8 8 128 -M REL 1e-3 -a \n");
    printf("	hpez -f -i test.dat    -o test.dat.qoz.out -3 8 8 128 -M ABS_AND_REL -A 1 -R 1e-3 -a -q 4\n");
    printf("	hpez -f -i test.dat    -o test.dat.qoz.out -3 8 8 128 -c qoz.config -q 2\n");
    printf("	hpez -f -i test.dat    -o test.dat.qoz.out -3 8 8 128 -c qoz.config -M ABS 1e-3 -a -q 1\n");
    exit(0);
}

void usage_sz2() {
    printf("Note: below are the supported command line arguments in SZ2 style\n");
    printf("Usage: qoz <options>\n");
    printf("Options:\n");
    printf("* operation type:\n");
    printf("	-z <compressed file>: the compression operation with an optionally specified output file.\n");
    printf("                          (the compressed file will be named as <input_file>.qoz if not specified)\n");
    printf("	-x <decompressed file>: the decompression operation with an optionally specified output file\n");
    printf("                      (the decompressed file will be named as <cmpred_file>.out if not specified)\n");
//    printf("	-p: print meta data (configuration info)\n");
    printf("	-h: print the help information\n");
    printf("	-v: print the version number\n");
    printf("* data type:\n");
    printf("	-f: single precision (float type)\n");
    printf("	-d: double precision (double type)\n");
    printf("* configuration file: \n");
    printf("	-c <configuration file> : configuration file qoz.config\n");
    printf("* error control: (the error control parameters here will overwrite the setting in qoz.config)\n");
    printf("	-M <error bound mode> : 10 options as follows. \n");
    printf("		ABS (absolute error bound)\n");
    printf("		REL (value range based error bound, so a.k.a., VR_REL)\n");
    printf("		ABS_AND_REL (using min{ABS, REL})\n");
    printf("		ABS_OR_REL (using max{ABS, REL})\n");
    printf("		PSNR (peak signal-to-noise ratio)\n");
    printf("		NORM (norm2 error : sqrt(sum(xi-xi')^2)\n");
//    printf("		PW_REL (point-wise relative error bound)\n");
    printf("	-A <absolute error bound>: specifying absolute error bound\n");
    printf("	-R <value_range based relative error bound>: specifying relative error bound\n");
//    printf("	-P <point-wise relative error bound>: specifying point-wise relative error bound\n");
    printf("	-S <PSNR>: specifying PSNR\n");
    printf("	-N <normErr>: specifying normErr\n");
    printf("* input data file:\n");
    printf("	-i <original data file> : original data file\n");
    printf("	-s <compressed data file> : compressed data file in decompression\n");
    printf("* output type of decompressed file: \n");
    printf("	-b (by default) : decompressed file stored in binary format\n");
    printf("	-t : decompreadded file stored in text format\n");
//    printf("	-T : pre-processing with Tucker Tensor Decomposition\n");
    printf("* dimensions: \n");
    printf("	-1 <nx> : dimension for 1D data such as data[nx]\n");
    printf("	-2 <nx> <ny> : dimensions for 2D data such as data[ny][nx]\n");
    printf("	-3 <nx> <ny> <nz> : dimensions for 3D data such as data[nz][ny][nx] \n");
    printf("	-4 <nx> <ny> <nz> <np>: dimensions for 4D data such as data[np][nz][ny][nx] \n");
    printf("* print compression results: \n");
    printf("	-a : print compression results such as distortions\n");
    printf("* examples: \n");
    printf("	qoz -z -f -c qoz.config -i testdata/x86/testfloat_8_8_128.dat -3 8 8 128\n");
    printf("	qoz -z -f -c qoz.config -M ABS -A 1E-3 -i testdata/x86/testfloat_8_8_128.dat -3 8 8 128\n");
    printf("	qoz -x -f -s testdata/x86/testfloat_8_8_128.dat.qoz -3 8 8 128\n");
    printf("	qoz -x -f -s testdata/x86/testfloat_8_8_128.dat.qoz -i testdata/x86/testfloat_8_8_128.dat -3 8 8 128 -a\n");
    printf("	qoz -z -d -c qoz.config -i testdata/x86/testdouble_8_8_128.dat -3 8 8 128\n");
    printf("	qoz -x -d -s testdata/x86/testdouble_8_8_128.dat.qoz -3 8 8 128\n");
    printf("	qoz -p -s testdata/x86/testdouble_8_8_128.dat.qoz\n");
    exit(0);
}

template<class T>
void compress(std::array<char *,3>inPaths, std::array<char *,3>cmpPaths, QoZ::Config &conf) {//conf changed to reference
    std::array<T*,3> data;
    for (auto i:{0,1,2}){
        data[i] = new T[conf.num];
        QoZ::readfile<T>(inPaths[i], conf.num, data[i]);
    }

    std::array<size_t,3> outSizes;
    QoZ::Timer timer(true);
    auto bytes = SZ_compress<T>(conf, data, outSizes);
    double compress_time = timer.stop();
   

    char outputFilePath[1024];
    for (auto i:{0,1,2}){
        if (cmpPaths[i] == nullptr) {
            sprintf(outputFilePath, "%s.hpez", inPaths[i]);
        } else {
            strcpy(outputFilePath, cmpPaths[i]);
        }
    
    
   
        QoZ::writefile(outputFilePath, bytes[i], outSizes[i]);


        
        
        printf("File %s:\n",inPaths[i]);
        printf("compression ratio = %.2f \n", conf.num * 1.0 * sizeof(T) / outSizes[i]);
        printf("compressed data file = %s\n", outputFilePath);
       

        
        delete[] data[i];
        delete[] bytes[i];
    }
    printf("compression time = %f\n", compress_time);

  

   
}

template<class T>
void decompress(std::array<char*,3> inPaths, std::array<char*,3> cmpPaths,  std::array<char*,3> decPaths,
                QoZ::Config &conf,
                int binaryOutput, int printCmpResults) {//conf changed to reference

    std::array<size_t,3> cmpSizes;

    std::array<char*,3> cmpData;

    for(auto i:{0,1,2}){
        cmpData[i]=QoZ::readfile<char>(cmpPaths[i], cmpSizes[i]).get();
        int confSize;
        std::cout<<cmpSizes[i]<<std::endl;
        memcpy(&confSize, cmpData[i] + (cmpSizes[i] - sizeof(int)), sizeof(int));
        std::cout<<confSize<<std::endl;
    }

    /*
    auto cmpData1 = QoZ::readfile<char>(cmpPath1, cmpSize1);
    auto cmpData2 = QoZ::readfile<char>(cmpPath2, cmpSize2); 
    auto cmpData3 = QoZ::readfile<char>(cmpPath3, cmpSize3); 
    */    

    QoZ::Timer timer(true);
    std::array<T *,3>decData = SZ_decompress<T>(conf, cmpData, cmpSizes);
    double decompress_time = timer.stop();

    char outputFilePaths[3][1024];
    for(auto i:{0,1,2}){
        if (decPaths[i] == nullptr) {
            sprintf(outputFilePaths[i], "%s.out", cmpPaths[i]);
        } else {
            strcpy(outputFilePaths[i], decPaths[i]);
        }
    }

    if (binaryOutput == 1) {
        for(auto i:{0,1,2})
            QoZ::writefile<T>(outputFilePaths[i], decData[i], conf.num);
        
    } else {
        for(auto i:{0,1,2})
            QoZ::writeTextFile<T>(outputFilePaths[i], decData[i], conf.num);
       
    }
    if (printCmpResults) {
        //compute the distortion / compression errors...
        size_t totalNbEle;
        std::array<T*,3> ori_data;
        for(auto i:{0,1,2}){
            ori_data[i] = QoZ::readfile<T>(inPaths[i], totalNbEle).get();
            
            assert(totalNbEle == conf.num);
        }
        //QoZ::verify<T>(ori_data.get(), decData, conf.num);
        //QoZ::verifyQoI<T>(ori_data.get(), decData, conf.dims, conf.qoiRegionSize);



        QoZ::verifyQoI_new<T>(ori_data, decData, conf);
    }
    for(auto i:{0,1,2}){
        delete[]   decData[i];
        printf("Compressed file %s:\n",cmpPaths[i]);
        printf("compression ratio = %f\n", conf.num * sizeof(T) * 1.0 / cmpSizes[i]);
        printf("decompressed file = %s\n", outputFilePaths[i]);
    

   
    }
    printf("decompression time = %f seconds.\n", decompress_time);
}

int main(int argc, char *argv[]) {
    bool binaryOutput = true;
    int printCmpResults = 0;
    bool compression = false;
    bool decompression = false;
    int dataType = SZ_FLOAT;
    char * inPath1 = nullptr;
    char * inPath2 = nullptr;
    char * inPath3 = nullptr;
    char * cmpPath1 = nullptr;
    char * cmpPath2 = nullptr;
    char * cmpPath3 = nullptr;
    char *conPath = nullptr;
    char *decPath1 = nullptr;
    char *decPath2 = nullptr;
    char *decPath3 = nullptr;
    bool delCmpPath = false;

    char *errBoundMode = nullptr;
    char *errBound = nullptr;
    char *absErrorBound = nullptr;
    char *relErrorBound = nullptr;
    char *pwrErrorBound = nullptr;
    char *psnrErrorBound = nullptr;
    char *normErrorBound = nullptr;
    char *tuningTarget = nullptr;
    int maxStep=0;
    int sampleBlockSize=0;

    bool sz2mode = false;
    int qoz=-1;
    bool testLorenzo=false;

    size_t r4 = 0;
    size_t r3 = 0;
    size_t r2 = 0;
    size_t r1 = 0;

    size_t i = 0;
    int status;
    if (argc == 1)
        usage();
    int width = -1;

    double quantile = -1;

    for (i = 1; i < argc; i++) {
        if (argv[i][0] != '-' || argv[i][2]) {
            if (argv[i][1] == 'h' && argv[i][2] == '2') {
                usage_sz2();
            } else {
                usage();
            }
        }
        switch (argv[i][1]) {
            case 'h':
                usage();
                exit(0);
            case 'v':
                printf("version: %s\n", QoZ_VER);
                exit(0);
            case 'b':
                binaryOutput = true;
                break;
            case 't':
                binaryOutput = false;
                break;
            case 'a':
                printCmpResults = 1;
                break;
            case 'z':
                compression = true;


                if (i + 1 < argc) {
                    cmpPath1 = argv[i + 1];
                    if (cmpPath1[0] != '-')
                        i++;
                    else
                        cmpPath1 = nullptr;
                }
                if (i + 1 < argc) {
                    cmpPath2 = argv[i + 1];
                    if (cmpPath2[0] != '-')
                        i++;
                    else
                        cmpPath2 = nullptr;
                }
                if (i + 1 < argc) {
                    cmpPath3 = argv[i + 1];
                    if (cmpPath3[0] != '-')
                        i++;
                    else
                        cmpPath3 = nullptr;
                }
                break;
            case 'x':
                sz2mode = true;
                decompression = true;
                if (i + 1 < argc) {
                    decPath1 = argv[i + 1];
                    if (decPath1[0] != '-')
                        i++;
                    else
                        decPath1 = nullptr;
                }
                if (i + 1 < argc) {
                    decPath2 = argv[i + 1];
                    if (decPath2[0] != '-')
                        i++;
                    else
                        decPath2 = nullptr;
                }
                if (i + 1 < argc) {
                    decPath3 = argv[i + 1];
                    if (decPath3[0] != '-')
                        i++;
                    else
                        decPath3 = nullptr;
                }
                break;
            case 'f':
                dataType = SZ_FLOAT;
                break;
            case 'd':
                dataType = SZ_DOUBLE;
                break;

            case 'I':
                if (++i == argc || sscanf(argv[i], "%d", &width) != 1) {
                    usage();
                }
                if (width == 32) {
                    dataType = SZ_INT32;
                } else if (width == 64) {
                    dataType = SZ_INT64;
                } else {
                    usage();
                }
                break;
            case 'i':
                if (i+3 >= argc || argv[i+1][0]== '-' || argv[i+2][0]=='-'|| argv[i+3][0]=='-')
                    usage();
                inPath1 = argv[i+1];
                inPath2 = argv[i+2];
                inPath3 = argv[i+3];
                i+=3;
                break;
            case 'q':
                if (++i == argc || sscanf(argv[i], "%d", &qoz) != 1)
                    usage();
                break;
            case 'l':
                testLorenzo = true;
                break;
            case 'o':
                if (i+3 >= argc || argv[i+1][0] == '-'|| argv[i+2][0] == '-'|| argv[i+3][0] == '-')
                    usage();
                decPath1 = argv[i+1];
                decPath2 = argv[i+2];
                decPath3 = argv[i+3];
                i+=3;
                break;
            case 's':
                sz2mode = true;
                if (i+3 >= argc|| argv[i+1][0]== '-' || argv[i+2][0]=='-'|| argv[i+3][0]=='-')
                    usage();

                cmpPath1 = argv[i + 1];
                cmpPath2 = argv[i + 2];
                cmpPath3 = argv[i + 3];
                i+=3;
                break;
            case 'c':
                if (++i == argc)
                    usage();
                conPath = argv[i];
                break;
            case '1':
                if (++i == argc || sscanf(argv[i], "%zu", &r1) != 1)
                    usage();
                break;
            case '2':
                if (++i == argc || sscanf(argv[i], "%zu", &r1) != 1 ||
                    ++i == argc || sscanf(argv[i], "%zu", &r2) != 1)
                    usage();
                break;
            case '3':
                if (++i == argc || sscanf(argv[i], "%zu", &r1) != 1 ||
                    ++i == argc || sscanf(argv[i], "%zu", &r2) != 1 ||
                    ++i == argc || sscanf(argv[i], "%zu", &r3) != 1)
                    usage();
                break;
            case '4':
                if (++i == argc || sscanf(argv[i], "%zu", &r1) != 1 ||
                    ++i == argc || sscanf(argv[i], "%zu", &r2) != 1 ||
                    ++i == argc || sscanf(argv[i], "%zu", &r3) != 1 ||
                    ++i == argc || sscanf(argv[i], "%zu", &r4) != 1)
                    usage();
                break;
            case 'M':
                if (++i == argc)
                    usage();
                errBoundMode = argv[i];
                if (i + 1 < argc && argv[i + 1][0] != '-') {
                    errBound = argv[++i];
                }
                break;
            case 'A':
                if (++i == argc)
                    usage();
                absErrorBound = argv[i];
                break;
            case 'R':
                if (++i == argc)
                    usage();
                relErrorBound = argv[i];
                break;
//            case 'P':
//                if (++i == argc)
//                    usage();
//                pwrErrorBound = argv[i];
//                break;
            case 'N':
                if (++i == argc)
                    usage();
                normErrorBound = argv[i];
                break;
            case 'S':
                if (++i == argc)
                    usage();
                psnrErrorBound = argv[i];
                break;

            case 'T':
                if (++i == argc)
                    usage();
                tuningTarget = argv[i];
                break;
            case 'C':
                if (++i == argc || sscanf(argv[i], "%d", &maxStep) != 1)
                        usage();
                break;
            case 'B':
                if (++i == argc || sscanf(argv[i], "%d", &sampleBlockSize) != 1)
                        usage();
                break;

            case 'u':
                if (++i == argc)
                    usage();
                quantile = atof (argv[i]);
                break;


            default:
                usage();
                break;
        }
    }

    if ((inPath3 == nullptr) && (cmpPath3 == nullptr)) {
        printf("Error: you need to specify either a raw binary data file or a compressed data file as input\n");
        usage();
        exit(0);
    }

    if (!sz2mode && inPath3 != nullptr && cmpPath3 != nullptr) {
        compression = true;
    }
    if (cmpPath3 != nullptr && decPath3 != nullptr) {
        decompression = true;
    }
    char cmpPathTmp1[1024],cmpPathTmp2[1024],cmpPathTmp3 [1024];
    if (inPath3 != nullptr && cmpPath3 == nullptr && decPath3 != nullptr) {
        compression = true;
        decompression = true;
        sprintf(cmpPathTmp1, "%s.hpez.tmp", inPath1);
        sprintf(cmpPathTmp2, "%s.hpez.tmp", inPath2);
        sprintf(cmpPathTmp3, "%s.hpez.tmp", inPath3);
        cmpPath1 = cmpPathTmp1;
        cmpPath2 = cmpPathTmp2;
        cmpPath3 = cmpPathTmp3;
        delCmpPath = true;
    }
    if (inPath3 == nullptr||errBoundMode == nullptr) {
        compression = false;
    }
    if (!compression && !decompression) {
        usage();
        exit(0);
    }

    QoZ::Config conf;
    if (r2 == 0) {
        conf = QoZ::Config(r1);
    } else if (r3 == 0) {
        conf = QoZ::Config(r2, r1);
    } else if (r4 == 0) {
        conf = QoZ::Config(r3, r2, r1);
    } else {
        conf = QoZ::Config(r4, r3, r2, r1);
    }
    if (qoz>=0){
        conf.QoZ=qoz;
    }
    if (testLorenzo){
        conf.testLorenzo=1;
    }
    if(maxStep>0)
        conf.maxStep=maxStep;
    if(conf.sampleBlockSize>0)
        conf.sampleBlockSize=sampleBlockSize;
    if (compression && conPath != nullptr) {
        conf.loadcfg(conPath);
    }


    if (errBoundMode != nullptr) {
        {
            // backward compatible with SZ2
            if (relErrorBound != nullptr) {
                conf.relErrorBound = atof(relErrorBound);
            }
            if (absErrorBound != nullptr) {
                conf.absErrorBound = atof(absErrorBound);
            }
            if (psnrErrorBound != nullptr) {
                conf.psnrErrorBound = atof(psnrErrorBound);
            }
            if (normErrorBound != nullptr) {
                conf.l2normErrorBound = atof(normErrorBound);
            }
        }
        if (strcmp(errBoundMode, QoZ::EB_STR[QoZ::EB_ABS]) == 0) {
            conf.errorBoundMode = QoZ::EB_ABS;
            if (errBound != nullptr) {
                conf.absErrorBound = atof(errBound);
            }
        } else if (strcmp(errBoundMode, QoZ::EB_STR[QoZ::EB_REL]) == 0 || strcmp(errBoundMode, "VR_REL") == 0) {
            conf.errorBoundMode = QoZ::EB_REL;
            if (errBound != nullptr) {
                conf.relErrorBound = atof(errBound);
            }
        } else if (strcmp(errBoundMode, QoZ::EB_STR[QoZ::EB_PSNR]) == 0) {
            conf.errorBoundMode = QoZ::EB_PSNR;
            if (errBound != nullptr) {
                conf.psnrErrorBound = atof(errBound);
            }
        } else if (strcmp(errBoundMode, QoZ::EB_STR[QoZ::EB_L2NORM]) == 0) {
            conf.errorBoundMode = QoZ::EB_L2NORM;
            if (errBound != nullptr) {
                conf.l2normErrorBound = atof(errBound);
            }
        } else if (strcmp(errBoundMode, QoZ::EB_STR[QoZ::EB_ABS_AND_REL]) == 0) {
            conf.errorBoundMode = QoZ::EB_ABS_AND_REL;
        } else if (strcmp(errBoundMode, QoZ::EB_STR[QoZ::EB_ABS_OR_REL]) == 0) {
            conf.errorBoundMode = QoZ::EB_ABS_OR_REL;
        } else {
            printf("Error: wrong error bound mode setting by using the option '-M'\n");
            usage();
            exit(0);
        }
    }

    if (tuningTarget!= nullptr) {
       
        if (strcmp(tuningTarget, "PSNR") == 0) {
            conf.tuningTarget = QoZ::TUNING_TARGET_RD;
        }
        else if (strcmp(tuningTarget, "CR") == 0) {
            conf.tuningTarget = QoZ::TUNING_TARGET_CR;
        }
        else if (strcmp(tuningTarget, "SSIM") == 0) {
            conf.tuningTarget = QoZ::TUNING_TARGET_SSIM;
        }
        else if (strcmp(tuningTarget, "AC") == 0) {
            conf.tuningTarget = QoZ::TUNING_TARGET_AC;
        }
        else {
            printf("Error: wrong tuning target setting by using the option '-T'\n");
            usage();
            exit(0);
        }
        
    }

    if (quantile>=0)
        conf.quantile = quantile;



    if (compression) {

        if (dataType == SZ_FLOAT) {
            compress<float>(std::array<char*,3>{inPath1, inPath2, inPath3}, std::array<char*,3>{cmpPath1, cmpPath2, cmpPath3}, conf);
        } 
        /*else if (dataType == SZ_DOUBLE) {
            compress<double>(inPath, cmpPath, conf);
        
        } 
        
        else if (dataType == SZ_INT32) {
            compress<int32_t>(inPath, cmpPath, conf);
        } else if (dataType == SZ_INT64) {
            compress<int64_t>(inPath, cmpPath, conf);
        }
        */
       
         else {
            printf("Error: data type not supported \n");
            usage();
            exit(0);
        }
        
    }
  
    if (decompression) {
        if (printCmpResults && inPath3 == nullptr) {
            printf("Error: Since you add -a option (analysis), please specify the original data path by -i <path>.\n");
            exit(0);
        }
       
        if (dataType == SZ_FLOAT) {
            decompress<float>(std::array<char*,3>{inPath1,inPath2,inPath3}, std::array<char*,3>{cmpPath1,cmpPath2,cmpPath3}, std::array<char*,3>{decPath1,decPath2,decPath3}, conf, binaryOutput, printCmpResults);

        } 
        /*else if (dataType == SZ_DOUBLE) {
            decompress<double>(inPath, cmpPath, decPath, conf, binaryOutput, printCmpResults);

        } 
        
        else if (dataType == SZ_INT32) {
            decompress<int32_t>(inPath, cmpPath, decPath, conf, binaryOutput, printCmpResults);
        } else if (dataType == SZ_INT64) {
            decompress<int64_t>(inPath, cmpPath, decPath, conf, binaryOutput, printCmpResults);
        }
         */
        else {
            printf("Error: data type not supported \n");
            usage();
            exit(0);
        }
    }

    if (delCmpPath) {
        remove(cmpPath1);
        remove(cmpPath2);
        remove(cmpPath3);
    }
    return 0;
}
