#ifndef _SZ_RUNLENGTH_ENCODER_HPP
#define _SZ_RUNLENGTH_ENCODER_HPP

#include "Encoder.hpp"
#include "QoZ/utils/MemoryUtil.hpp"
#include "QoZ/def.hpp"
#include <vector>

namespace QoZ {

    template<class T>
    class RunlengthEncoder : public concepts::EncoderInterface<T> {
    public:

        ~RunlengthEncoder() = default;

        void preprocess_encode(const std::vector<T> &bins, int stateNum) {
        };

        void preprocess_encode() {
        };

        size_t encode(const T * bins, size_t num_elements, uchar *&bytes) {
            int max = 0;
            size_t s = 0;
            for (size_t i = 1; i < num_elements; i++) {
                if (bins[i] != bins[i - 1]) {
                    write(bins[i - 1], bytes);
                    write(int(i - s), bytes);
                    if (int(i - s) > max) {
                        max = int(i - s);
                    }
                    s = i;
                }
            }
            write(bins[num_elements - 1], bytes);
            write(int(num_elements - s), bytes);
            printf("RunLengthEncoder max length = %d\n", max);
            return 0;
        };

        size_t encode(const std::vector<T> &bins, size_t num_elements, uchar *&bytes) {
            return encode(bins.data(), num_elements, bytes);
        }

        size_t encode(const std::vector<T> &bins, uchar *&bytes) {
            return encode(bins.data(), bins.size(), bytes);
        }

        void postprocess_encode() {};

        void preprocess_decode() {};

        std::vector<T> decode(const uchar *&bytes, size_t targetLength) {
            std::vector<T> bins(targetLength, 0);
            T value;
            int cnt;
            for (size_t i = 0; i < bins.size();) {
                read(value, bytes);
                read(cnt, bytes);
                for (size_t j = i; j < i + cnt; j++) {
                    bins[j] = value;
                }
                i += cnt;
            }
            return bins;
        };

        void postprocess_decode() {};

        uint save(uchar *&c) {
            return 0;
        };

        size_t size_est() {
            return 0;
        }

        void load(const uchar *&c, size_t &remaining_length) {};

    };
}
#endif
