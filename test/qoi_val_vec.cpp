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
    printf("Usage: qoi_val <options>\n");
    printf("Options:\n");
    printf("* input and output:\n");
    printf("    -i <path> : original 3 binary input files\n");
    printf("    -o <path> : 3 compressed output files, default in binary format\n");
//    printf("  -p: print meta data (configuration info)\n");
    printf("* data type:\n");
    printf("    -f: single precision (float type)\n");
    printf("    -d: double precision (double type)\n");
    printf("    -I <width>: integer type (width = 32 or 64)\n");
    printf("* configuration file: \n");
    printf("    -c <configuration file> : configuration file qoz.config\n");
    printf("* dimensions: \n");
    printf("    -1 <nx> : dimension for 1D data such as data[nx]\n");
    printf("    -2 <nx> <ny> : dimensions for 2D data such as data[ny][nx]\n");
    printf("    -3 <nx> <ny> <nz> : dimensions for 3D data such as data[nz][ny][nx] \n");
    printf("    -4 <nx> <ny> <nz> <np>: dimensions for 4D data such as data[np][nz][ny][nx] \n");
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
void qoiValidation(std::array<char*,3> inPaths, std::array<char*,3> decPaths,
                QoZ::Config &conf) {

        size_t totalNbEle;
        std::array<T*,3> ori_data;
        std::array<T*,3> dec_data;
        auto o1 = QoZ::readfile<T>(inPaths[0], totalNbEle);
        auto o2 = QoZ::readfile<T>(inPaths[1], totalNbEle);
        auto o3 = QoZ::readfile<T>(inPaths[2], totalNbEle);
       
        ori_data[0] = o1.get();
        ori_data[1] = o2.get();
        ori_data[2] = o3.get();

        auto d1 = QoZ::readfile<T>(decPaths[0], totalNbEle);
        auto d2 = QoZ::readfile<T>(decPaths[1], totalNbEle);
        auto d3 = QoZ::readfile<T>(decPaths[2], totalNbEle);
       
        dec_data[0] = d1.get();
        dec_data[1] = d2.get();
        dec_data[2] = d3.get();
            
        assert(totalNbEle == conf.num);
        
        //QoZ::verify<T>(ori_data.get(), decData, conf.num);
        //QoZ::verifyQoI<T>(ori_data.get(), decData, conf.dims, conf.qoiRegionSize);


        //conf.qoi = 1;//Forced verify qoi
        QoZ::verifyQoI_new<T>(ori_data, decData, conf);
    
}

int main(int argc, char *argv[]) {
   
    int dataType = SZ_FLOAT;
    char * inPath1 = nullptr;
    char * inPath2 = nullptr;
    char * inPath3 = nullptr;
    
    char *conPath = nullptr;
    char *decPath1 = nullptr;
    char *decPath2 = nullptr;
    char *decPath3 = nullptr;

    size_t r4 = 0;
    size_t r3 = 0;
    size_t r2 = 0;
    size_t r1 = 0;

    size_t i = 0;
    int status;
    if (argc == 1)
        usage();
    int width = -1;

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
           
                printCmpResults = 1;
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
            
            case 'o':
                if (i+3 >= argc || argv[i+1][0] == '-'|| argv[i+2][0] == '-'|| argv[i+3][0] == '-')
                    usage();
                decPath1 = argv[i+1];
                decPath2 = argv[i+2];
                decPath3 = argv[i+3];
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

            default:
                usage();
                break;
        }
    }

    if ((inPath1 == nullptr) ||(inPath2 == nullptr) ||(inPath3 == nullptr) || (decPath1 == nullptr)|| (decPath2 == nullptr)|| (decPath3 == nullptr)) {
        printf("Error: ori or dec data missing\n");
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
    
    if (compression && conPath != nullptr) {
        conf.loadcfg(conPath);
    }

  

    if (dataType == SZ_FLOAT) {
        qoiValidation<float>(std::array<char*,3>{inPath1, inPath2, inPath3}, std::array<char*,3>{decPath1, decPath2, decPath3}, conf);
    } 
    else if (dataType == SZ_DOUBLE) {
        qoiValidation<double>(std::array<char*,3>{inPath1, inPath2, inPath3}, std::array<char*,3>{decPath1, decPath2, decPath3}, conf);
        
    } 
        
    else if (dataType == SZ_INT32) {
        qoiValidation<int32_t>(std::array<char*,3>{inPath1, inPath2, inPath3}, std::array<char*,3>{decPath1, decPath2, decPath3}, conf);
    } else if (dataType == SZ_INT64) {
        qoiValidation<int64_t>(std::array<char*,3>{inPath1, inPath2, inPath3}, std::array<char*,3>{decPath1, decPath2, decPath3}, conf);
    }
    return 0;
}
