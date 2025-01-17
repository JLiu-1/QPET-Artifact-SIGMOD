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
    printf("* input:\n");
    printf("	-i <path> : original binary input file\n");
    printf("	-o <path> : decompressed output file, default in binary format\n");
//    printf("	-p: print meta data (configuration info)\n");
    printf("* data type:\n");
    printf("	-f: single precision (float type)\n");
    printf("	-d: double precision (double type)\n");
    printf("	-I <width>: integer type (width = 32 or 64)\n");
    printf("* configuration file: \n");
    printf("	-c <configuration file> : configuration file hpez.config\n");
    printf("* dimensions: \n");
    printf("	-1 <nx> : dimension for 1D data such as data[nx]\n");
    printf("	-2 <nx> <ny> : dimensions for 2D data such as data[ny][nx]\n");
    printf("	-3 <nx> <ny> <nz> : dimensions for 3D data such as data[nz][ny][nx] \n");
    printf("	-4 <nx> <ny> <nz> <np>: dimensions for 4D data such as data[np][nz][ny][nx] \n");
    printf("* output: Validation results of decompression data/Qoi quality\n");
    printf("Example: qoi_val -f -3 128 128 128 -c qoi.config -i ori.dat -o dec.dat");
    exit(0);
}


template<class T>
void qoiValidation(char *inPath, char *decPath,
                QoZ::Config &conf) {

   
    //compute the distortion / compression errors...
    size_t totalNbEle;
    std::cout<<"a1"<<std::endl;
    auto ori_data = QoZ::readfile<T>(inPath, totalNbEle);
     std::cout<<"a2"<<std::endl;
    assert(totalNbEle == conf.num);
    auto dec_data = QoZ::readfile<T>(decPath, totalNbEle);
     std::cout<<"a3"<<std::endl;
    //QoZ::verify<T>(ori_data.get(), decData, conf.num);
    //QoZ::verifyQoI<T>(ori_data.get(), decData, conf.dims, conf.qoiRegionSize);


   
    QoZ::verifyQoI_new<T>(ori_data.get(), dec_data.get(), conf);
     std::cout<<"a4"<<std::endl;
    

}

int main(int argc, char *argv[]) {

    int dataType = SZ_FLOAT;
    char *inPath = nullptr;
    char *conPath = nullptr;
    char *decPath = nullptr;
    

  

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
            
            usage();
            
        }
        switch (argv[i][1]) {
            case 'h':
                usage();
                exit(0);
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
                if (++i == argc)
                    usage();
                inPath = argv[i];
                break;
            case 'o':
                if (++i == argc)
                    usage();
                decPath = argv[i];
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

    if ((inPath == nullptr) || (decPath == nullptr)) {
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
    if (conPath != nullptr) {
        conf.loadcfg(conPath);
    }

    if (dataType == SZ_FLOAT) {
        qoiValidation<float>(inPath, decPath, conf);
    } 
    else if (dataType == SZ_DOUBLE) {
        qoiValidation<double>(inPath, decPath, conf);
        
    } 
        
    else if (dataType == SZ_INT32) {
        qoiValidation<int32_t>(inPath, decPath, conf);
    } else if (dataType == SZ_INT64) {
        qoiValidation<int64_t>(inPath, decPath, conf);
    }


    

    return 0;
}
