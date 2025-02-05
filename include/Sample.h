

#ifndef SAMPLE_H
#define SAMPLE_H

namespace sperr {

    /*
    template<class T, uint N>
    inline void
    sample_block_4d(T *data, std::vector<T> & sampling_data , std::vector<size_t> &dims, std::vector<size_t> &starts,size_t block_size) {
        assert(dims.size() == N);
        assert(starts.size() == N);
        
        
        size_t sample_num = block_size*block_size*block_size;
        //std::vector<T> sampling_data(sample_num, 0);

       
//        auto sampling_time = timer.stop();
//        printf("Generate sampling data, block = %lu percent = %.3f%% Time = %.3f \n", sampling_block, sample_num * 100.0 / num,
//               sampling_time);
        //return sampling_data;
    }
 
    template<class T, uint N>
    inline void
    calculate_interp_error_vars(T *data, std::vector<size_t> &dims,std::vector<double> &vars,uint8_t interp_op=0,uint8_t nat=0, size_t stride=8,size_t interp_stride=1,T abs_eb=0.0){

        vars.resize(N,0);
        size_t count=0;
        if(stride<2)
            stride=2;
        if(N==3){
            size_t dimx=dims[0],dimy=dims[1],dimz=dims[2],dimyz=dimy*dimz;
            size_t is1=dimyz*interp_stride,is3x1=3*is1,is2=dimz*interp_stride,is3x2=3*is2,is3=interp_stride,is3x3=3*is3;
            for (size_t i = 3*interp_stride; i+3*interp_stride < dimx; i+=(stride/2)*2*interp_stride) {
                for (size_t j = 3*interp_stride; j+3*interp_stride < dimy; j+=(stride/2)*2*interp_stride) {
                    for (size_t k = 3*interp_stride; k+3*interp_stride < dimz; k+=(stride/2)*2*interp_stride) {
                        count+=1;
                        size_t idx=i*dimyz+j*dimz+k;
                        T *d= data+idx;
                        T cur_value=*d;
                        if(interp_op==1){
                            auto interp_cubic=nat?interp_cubic_2<T>:interp_cubic_1<T>;
                            T interp_err=interp_cubic(*(d - is3x1), *(d - is1), *(d + is1), *(d + is3x1))-cur_value;
                            vars[0]+=interp_err*interp_err;
                            interp_err=interp_cubic(*(d - is3x2), *(d - is2), *(d +is2), *(d + is3x2))-cur_value;
                            vars[1]+=interp_err*interp_err;
                            interp_err=interp_cubic(*(d - is3x3), *(d - is3), *(d + is3), *(d + is3x3))-cur_value;
                            vars[2]+=interp_err*interp_err;
                        }
                        else{
                            T interp_err=interp_linear<T>( *(d -  is1), *(d + is1))-cur_value;
                            vars[0]+=interp_err*interp_err;
                            interp_err=interp_linear<T>( *(d -  is2), *(d + is2))-cur_value;
                            vars[1]+=interp_err*interp_err;
                            interp_err=interp_linear<T>( *(d -  is3), *(d + is3) )-cur_value;
                            vars[2]+=interp_err*interp_err;

                        }
                    }
                }
            }
        }
        else if(N==2){
            size_t  dimx=dims[0],dimy=dims[1];
            size_t is1=dimy*interp_stride,is3x1=3*is1,is2=interp_stride,is3x2=3*is2;
            for (size_t i = 3*interp_stride;i+3*interp_stride < dimx; i+=(stride/2)*2*interp_stride) {
                for (size_t j = 3*interp_stride; j+3*interp_stride < dimy; j+=(stride/2)*2*interp_stride) {
                 
                    count+=1;
                    size_t idx=i*dimy+j;
                    T *d= data+idx;
                    T cur_value=*d;
                    if(interp_op==1){
                        auto interp_cubic=nat?interp_cubic_2<T>:interp_cubic_1<T>;
                        T interp_err=interp_cubic(*(d - is3x1), *(d - is1), *(d + is1), *(d + is3x1))-cur_value;
                        vars[0]+=interp_err*interp_err;
                        interp_err=interp_cubic(*(d - is3x2), *(d - is2), *(d +is2), *(d + is3x2))-cur_value;
                        vars[1]+=interp_err*interp_err;
                    }
                    else{
                        T interp_err=interp_linear<T>( *(d -  is1), *(d + is1))-cur_value;
                        vars[0]+=interp_err*interp_err;
                        interp_err=interp_linear<T>( *(d -  is2), *(d + is2))-cur_value;
                        vars[1]+=interp_err*interp_err;
                    }
                }
            }

        }
        
        for (size_t i=0;i<N;i++){
            if(count>0)
                vars[i]/=double(count);
            else
                vars[i]=1.0;
            if(interp_op==1){
                if(nat)
                    vars[i]+=abs_eb*abs_eb*(1.0/12)*0.6725;
                else
                    vars[i]+=abs_eb*abs_eb*(1.0/12)*0.640625;
            }
            else{
                vars[i]+=abs_eb*abs_eb*(1.0/12)*0.5;
            }
        }
    }

    template<uint N>
    inline void
    preprocess_vars(std::vector<double>&vars){
        if(N==2){
            double a=vars[1],b=vars[0];
            vars[0]=a/(a+b);
            vars[1]=b/(a+b);
        }
        else if (N==3){
            double a=vars[1]*vars[2],b=vars[0]*vars[2],c=vars[0]*vars[1];
            vars[0]=a/(a+b+c);
            vars[1]=b/(a+b+c);
            vars[2]=c/(a+b+c);
        }
    }




template<class T, uint N>
    inline int
    calculate_interp_coeffs(T *data, std::vector<size_t> &dims,std::vector<double> &coeffs, size_t stride=2){

        std::vector<double>xs,ys; 
        size_t stride2x=2*stride;
        if(N==3){
            size_t dimx=dims[0],dimy=dims[1],dimz=dims[2],dimyz=dimy*dimz;
            for (size_t i = 3; i < dimx-3; i+=stride) {
                for (size_t j = 3; j < dimy-3; j+=stride) {
                    for (size_t k = 3; k < dimz-3; k+=stride) {
                        if(i%stride2x==0 and j%stride2x==0 and k%stride2x==0)
                            continue;
                        size_t idx=i*dimyz+j*dimz+k;
                        T *d= data+idx;
                        T cur_value=*d;
                        std::vector<double>temp_xs={*(d - 3*dimyz),*(d - dimyz),*(d + dimyz), *(d + 3*dimyz),*(d - 3*dimz),*(d - dimz),*(d + dimz),*(d + 3*dimz),*(d - 3), *(d - 1), *(d + 1), *(d + 3)};
                        std::vector<double>temp_ys={cur_value,cur_value,cur_value};
                        xs.insert(xs.end(),temp_xs.begin(),temp_xs.end());
                        ys.insert(ys.end(),temp_ys.begin(),temp_ys.end());                       

                    }
                }
            }
        }

        else if(N==2){
            size_t  dimx=dims[0],dimy=dims[1];
            for (size_t i = 3; i < dimx-3; i+=stride) {
                for (size_t j = 3; j < dimy-3; j+=stride) {
                    if(i%stride2x==0 and j%stride2x==0)
                        continue;
                 
                    //count+=1;
                    size_t idx=i*dimy+j;
                    T *d= data+idx;
                    T cur_value=*d;

                    std::vector<double>temp_xs={*(d - 3*dimy),*(d - dimy),*(d + dimy),*(d + 3*dimy),*(d - 3), *(d - 1), *(d + 1), *(d + 3)};
                    std::vector<double>temp_ys={cur_value,cur_value,cur_value};
                    xs.insert(xs.end(),temp_xs.begin(),temp_xs.end());
                    ys.insert(xs.end(),temp_ys.begin(),temp_ys.end());
                }
            }

        }
        int status;
        auto reg_res=QoZ::Regression(xs.data(),ys.size(),4,ys.data(),status);
        if(status==0){
            coeffs.resize(4);
            for(size_t i=0;i<4;i++)
                coeffs[i]=reg_res[i];
        }
        return status;
    }
    


*/


    template <class T>
    void blockwise_profiling(const T *data, const std::vector<size_t> &dims, const std::vector<size_t> &starts,const size_t &blocksize,double & mean,double & sigma2,double & range){
        size_t N=dims.size();

        if(N==2){
            size_t dimx=dims[0],dimy=dims[1],startx=starts[0],starty=starts[1],element_num=blocksize*blocksize;
            double sum=0,max,min;
            sigma2=0;
            size_t start_idx=startx*dimy+starty;
            max=min=data[start_idx];
            for(size_t i=startx;i<startx+blocksize;i++){
                for(size_t j=starty;j<starty+blocksize;j++){
                    size_t cur_idx=i*dimy+j;
                    T value=data[cur_idx];
                    sum+=value;
                    
                    max=value>max?value:max;
                    min=value<min?value:min;

                }

            }
            mean=sum/element_num;

            for(size_t i=startx;i<startx+blocksize;i++){
                for(size_t j=starty;j<starty+blocksize;j++){
                    size_t cur_idx=i*dimy+j;
                    T value=data[cur_idx];
                    sigma2+=(value-mean)*(value-mean);

                }

            }
            
            range=max-min;
            sigma2/=element_num;
        }


        else if(N==3){
            size_t dimx=dims[0],dimy=dims[1],dimz=dims[2],startx=starts[0],starty=starts[1],startz=starts[2],element_num=blocksize*blocksize*blocksize;
            size_t dimyz=dimy*dimz;
            double sum=0,max,min;
            sigma2=0;
            size_t start_idx=startx*dimyz+starty*dimz+startz;
            max=min=data[start_idx];
            for(size_t i=startx;i<startx+blocksize;i++){
                for(size_t j=starty;j<starty+blocksize;j++){
                    for(size_t k=startz;k<startz+blocksize;k++){
                        size_t cur_idx=i*dimyz+j*dimz+k;
                        T value=data[cur_idx];
                        sum+=value;
                       
                        max=value>max?value:max;
                        min=value<min?value:min;
                    }

                }

            }
            mean=sum/element_num;
            for(size_t i=startx;i<startx+blocksize;i++){
                for(size_t j=starty;j<starty+blocksize;j++){
                    for(size_t k=startz;k<startz+blocksize;k++){
                        size_t cur_idx=i*dimyz+j*dimz+k;
                        T value=data[cur_idx];
                        sigma2+=(value-mean)*(value-mean);

                        
                    }

                }

            }

           
            range=max-min;
            sigma2/=element_num;
        }
    }
    


    template<class T, uint N>
    inline void
    profiling_block_3d(T *data, std::vector<size_t> &dims, std::vector< std::vector<size_t> > &starts,size_t block_size, double abseb,size_t stride=4) {
        assert(dims.size() == N);
        if (stride==0)
            stride=block_size;
        size_t dimx=dims[0],dimy=dims[1],dimz=dims[2],dimyz=dimy*dimz;
        
        for (size_t i = 0; i < dimx-block_size; i+=block_size) {
            for (size_t j = 0; j < dimy-block_size; j+=block_size) {
                for (size_t k = 0; k < dimz-block_size; k+=block_size) {
                    //std::cout<<i<<" "<<j<<" "<<k<<std::endl;
                    size_t start_idx=i*dimyz+j*dimz+k;
                    T min=data[start_idx];
                    T max=data[start_idx];
                    for (int ii=0;ii<=block_size;ii+=stride){
                        for(int jj=0;jj<=block_size;jj+=stride){
                            for (int kk=0;kk<=block_size;kk+=stride){
                                size_t cur_idx=start_idx+ii*dimyz+jj*dimz+kk;
                                T cur_value=data[cur_idx];
                                if (cur_value<min)
                                    min=cur_value;
                                else if (cur_value>max)
                                    max=cur_value;

                            }
                        }
                    }
                    if (max-min>abseb){
                        size_t a[3]={i,j,k};
                        starts.push_back(std::vector<size_t>(a,a+3));
                    }


                    
                    
                }
            }
        }
        //current has a problem. May return no blocks. Thinking how to better solve it.
//        auto sampling_time = timer.stop();
//        printf("Generate sampling data, block = %lu percent = %.3f%% Time = %.3f \n", sampling_block, sample_num * 100.0 / num,
//               sampling_time);
       // return sampling_data;
    }
    


    template<class T, uint N>
    inline void
    sample_blocks(T *data, std::vector<T> & sampling_data, std::vector<size_t> &dims, std::vector<size_t> &starts,size_t block_size) {
        assert(dims.size() == N);
        assert(starts.size() == N);
        if(N==3){
        
            size_t sample_num = block_size*block_size*block_size;
            sampling_data.resize(sample_num, 0);

            size_t startx=starts[0],starty=starts[1],startz=starts[2],dimx=dims[0],dimy=dims[1],dimz=dims[2];
            size_t square_block_size=block_size*block_size,dimyz=dimy*dimz;
            for (size_t i = 0; i < block_size; i++) {
                for (size_t j = 0; j < block_size; j++) {
                    for (size_t k = 0; k < block_size; k++) {
                        size_t sample_idx=i*square_block_size+j*block_size+k;
                        size_t idx=(i+startx)*dimyz+(j+starty)*dimz+k+startz;
                        sampling_data[sample_idx]=data[idx];
                        
                    }
                }
            }
        }
        else if (N==2){
            size_t sample_num = block_size*block_size;
            sampling_data.resize(sample_num, 0);
            size_t startx=starts[0],starty=starts[1],dimx=dims[0],dimy=dims[1];
            
            for (size_t i = 0; i < block_size; i++) {
                for (size_t j = 0; j < block_size; j++) {
                    
                    size_t sample_idx=i*block_size+j;
                    size_t idx=(i+startx)*dimy+(j+starty);
                    sampling_data[sample_idx]=data[idx];
                        
                    
                }
            }

        }
        else if(N==1){
            size_t sample_num = block_size;
            sampling_data.resize(sample_num, 0);

            size_t startx=starts[0],dimx=dims[0];
            
            for (size_t i = 0; i < block_size; i++) {
                
                    
                size_t sample_idx=i;
                size_t idx=(i+startx);
                sampling_data[sample_idx]=data[idx];
                        
                    
                
            }

        }
//        auto sampling_time = timer.stop();
//        printf("Generate sampling data, block = %lu percent = %.3f%% Time = %.3f \n", sampling_block, sample_num * 100.0 / num,
//               sampling_time);
       // return sampling_data;
    }
 
    template<class T, uint N>
    inline void
    profiling_block_2d(T *data, std::vector<size_t> &dims, std::vector< std::vector<size_t> > &starts,size_t block_size, double abseb,size_t stride=4) {
        assert(dims.size() == N);
        if (stride==0)
            stride=block_size;
        
        size_t dimx=dims[0],dimy=dims[1];
        
        for (size_t i = 0; i < dimx-block_size; i+=block_size) {
            for (size_t j = 0; j < dimy-block_size; j+=block_size) {
                
                size_t start_idx=i*dimy+j;
                T min=data[start_idx];
                T max=data[start_idx];
                for (int ii=0;ii<=block_size;ii+=stride){
                    for(int jj=0;jj<=block_size;jj+=stride){
                           
                        size_t cur_idx=start_idx+ii*dimy+jj;
                        T cur_value=data[cur_idx];
                        if (cur_value<min)
                            min=cur_value;
                        else if (cur_value>max)
                            max=cur_value;

                    }
                }
                    
                if (max-min>abseb){
                     size_t a[2]={i,j};
                    starts.push_back(std::vector<size_t>(a,a+2));
                }


                    
                    
            }
        }
        
//        auto sampling_time = timer.stop();
//        printf("Generate sampling data, block = %lu percent = %.3f%% Time = %.3f \n", sampling_block, sample_num * 100.0 / num,
//               sampling_time);
       // return sampling_data;
    }

    template<class T, int N>
    void sampleBlocks(T *data,std::vector<size_t> &dims, size_t sampleBlockSize,std::vector< std::vector<T> > & sampled_blocks,double sample_rate,int profiling ,std::vector<std::vector<size_t> > &starts,int var_first=0){
        for(int i=0;i<sampled_blocks.size();i++){
                    std::vector< T >().swap(sampled_blocks[i]);                
                }
                std::vector< std::vector<T> >().swap(sampled_blocks);
        for(int i=0;i<sampled_blocks.size();i++){
            std::vector< T >().swap(sampled_blocks[i]);                  
        }
        std::vector< std::vector<T> >().swap(sampled_blocks);                               
        size_t totalblock_num=1;
        for(int i=0;i<N;i++){                        
            totalblock_num*=(int)((dims[i]-1)/sampleBlockSize);
        }               
        size_t idx=0,block_idx=0;   
        if(profiling){
            size_t num_filtered_blocks=starts.size();    
            if(var_first==0){  
                size_t sample_stride=(size_t)(num_filtered_blocks/(totalblock_num*sample_rate));
                if(sample_stride<=0)
                    sample_stride=1;
                
                for(size_t i=0;i<num_filtered_blocks;i+=sample_stride){
                    std::vector<T> s_block;
                    sample_blocks<T,N>(data, s_block,dims, starts[i],sampleBlockSize+1);
                    sampled_blocks.push_back(s_block);
                    
                }
                
            }
            else{
                std::vector< std::pair<double,std::vector<size_t> > >block_heap;
                for(size_t i=0;i<num_filtered_blocks;i++){
                    double mean,sigma2,range;
                    blockwise_profiling<T>(data,dims, starts[i],sampleBlockSize+1, mean,sigma2,range);
                    block_heap.push_back(std::pair<double,std::vector<size_t> >(sigma2,starts[i]));
                    
                }
                std::make_heap(block_heap.begin(),block_heap.end());
              

                size_t sampled_block_num=totalblock_num*sample_rate;
                if(sampled_block_num>num_filtered_blocks)
                    sampled_block_num=num_filtered_blocks;
                if(sampled_block_num==0)
                    sampled_block_num=1;

                for(size_t i=0;i<sampled_block_num;i++){
                    std::vector<T> s_block;
                 
                    sample_blocks<T,N>(data, s_block,dims, block_heap.front().second,sampleBlockSize+1);
                  
                    sampled_blocks.push_back(s_block);
                    std::pop_heap(block_heap.begin(),block_heap.end());
                    block_heap.pop_back();
                   
                }
            }
        }               
        else{
            if(var_first==0){
                size_t sample_stride=(size_t)(1.0/sample_rate);
                if(sample_stride<=0)
                    sample_stride=1;
                if (N==2){                        
                    for (size_t x_start=0;x_start<dims[0]-sampleBlockSize;x_start+=sampleBlockSize){                           
                        for (size_t y_start=0;y_start<dims[1]-sampleBlockSize;y_start+=sampleBlockSize){
                            if (idx%sample_stride==0){
                                std::vector<size_t> starts{x_start,y_start};
                                std::vector<T> s_block;
                                sample_blocks<T,N>(data, s_block,dims, starts,sampleBlockSize+1);
                                sampled_blocks.push_back(s_block);
                            }
                            idx+=1;
                        }
                    }
                }
                else if (N==3){                  
                    for (size_t x_start=0;x_start<dims[0]-sampleBlockSize;x_start+=sampleBlockSize){                          
                        for (size_t y_start=0;y_start<dims[1]-sampleBlockSize;y_start+=sampleBlockSize){
                            for (size_t z_start=0;z_start<dims[2]-sampleBlockSize;z_start+=sampleBlockSize){
                                if (idx%sample_stride==0){
                                    std::vector<size_t> starts{x_start,y_start,z_start};
                                    std::vector<T> s_block;
                                    sample_blocks<T,N>(data, s_block,dims, starts,sampleBlockSize+1);
                                    sampled_blocks.push_back(s_block);
                                }
                                idx+=1;
                            }
                        }
                    }
                }
            }
            else{
                std::vector <std::vector<size_t> > blocks_starts;
                if (N==2){  
                    for (size_t x_start=0;x_start<dims[0]-sampleBlockSize;x_start+=sampleBlockSize){                           
                        for (size_t y_start=0;y_start<dims[1]-sampleBlockSize;y_start+=sampleBlockSize){
                           
                                blocks_starts.push_back(std::vector<size_t>{x_start,y_start});
                        }
                    }

                }
                else if (N==3){           
                    for (size_t x_start=0;x_start<dims[0]-sampleBlockSize;x_start+=sampleBlockSize){                          
                        for (size_t y_start=0;y_start<dims[1]-sampleBlockSize;y_start+=sampleBlockSize){
                            for (size_t z_start=0;z_start<dims[2]-sampleBlockSize;z_start+=sampleBlockSize){
                                blocks_starts.push_back(std::vector<size_t>{x_start,y_start,z_start});
                            }
                        }
                    }
                

                    std::vector< std::pair<double,std::vector<size_t> > >block_heap;
                    for(size_t i=0;i<totalblock_num;i++){
                        double mean,sigma2,range;
                        blockwise_profiling<T>(data,dims, blocks_starts[i],sampleBlockSize+1, mean,sigma2,range);
                        block_heap.push_back(std::pair<double,std::vector<size_t> >(sigma2,blocks_starts[i]));
                    }
                    std::make_heap(block_heap.begin(),block_heap.end());
                    size_t sampled_block_num=totalblock_num*sample_rate;
                    if(sampled_block_num==0)
                        sampled_block_num=1;
                    for(size_t i=0;i<sampled_block_num;i++){
                        std::vector<T> s_block;
                        sample_blocks<T,N>(data, s_block,dims, block_heap.front().second,sampleBlockSize+1);
                        sampled_blocks.push_back(s_block);
                        std::pop_heap(block_heap.begin(),block_heap.end());
                        block_heap.pop_back();
                    }

                }
            }
        }
    }


   
    


};


#endif
