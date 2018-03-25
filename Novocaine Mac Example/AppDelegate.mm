// Copyright (c) 2012 Alex Wiltschko
// 
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:
// 
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.


#import "AppDelegate.h"

@implementation AppDelegate

- (void)dealloc
{
    if (_ringBufferIn){
        delete _ringBufferIn;
    }
    if (_ringBufferOut){
        delete _ringBufferOut;
    }
}

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification
{
    
    self.audioManager = [Novocaine audioManager];
    
    int ringBufferLength = 64*512;
    self.ringBufferIn = new RingBuffer(ringBufferLength, 2);
    self.ringBufferOut = new RingBuffer(ringBufferLength, 2);
    self.ringBufferOverlap = new RingBuffer(4*512, 2);
    
    __weak AppDelegate * wself = self;
 
   
    __block int dumpHit = 0;
    __block int dumpMatch = 10;
    __block int dumpCount = 1;

    // doing FFT's that are twice as long as 2*numFrames == 2*512
    
    int N = 1024;
    int log2n = 10;     // 2^10 = 1024
    
    // from 1024 containing both channels at once to two arrays of 512 for L and R separate
    float* deInterleavedL;
    float* deInterleavedR;
    deInterleavedL = new float[N];
    deInterleavedR = new float[N];
    
    // numFrames is 512, numChannels is 2
        
    [self.audioManager setInputBlock:^(float *data,
                                       UInt32 numFrames,
                                       UInt32 numChannels) {
        
        // 1. stash all input samples into an INPUT ring buffer
        /*
        float volume = 0.75;
        vDSP_vsmul(data, 1, &volume, data, 1, numFrames*numChannels);
        */
        
        for(int i=0; i<numFrames; i++){
            deInterleavedL[i] = data[2*i];
            // deInterleavedR[i] = data[2*i + 1];
        }
        
        wself.ringBufferIn->AddNewFloatData(deInterleavedL, numFrames, 0);
        // wself.ringBufferIn->AddNewFloatData(deInterleavedR, numFrames, 1);
        
//         wself.ringBufferIn->AddNewInterleavedFloatData(data, numFrames, numChannels);
        
        
//        if(!dumpHit){
//            if(dumpCount==dumpMatch) {
//                dumpHit = 1;
//                printf("hit it\n");
//                writeOutput(deInterleavedL, numFrames, 44100);
//            }
//            else
//                dumpCount++;
//        }
        
    }];
  
    


    // don't think i need a window
    float* myWindow;
    myWindow = new float[N/2];
    vDSP_hann_window(myWindow, N/2, vDSP_HANN_DENORM);
    
    FFTSetup setup = vDSP_create_fftsetup(log2n, kFFTRadix2);
    
    // doing overlap save, need to grab 1024 time samples each go round
    float* bufferInL;
    bufferInL = (float*)malloc(N*sizeof(float));
    float* bufferInR;
    bufferInR = (float*)malloc(N*sizeof(float));
    
    float* bufferOutL;
    bufferOutL = (float*)malloc(N*sizeof(float));
    float* bufferOutR;
    bufferOutR = (float*)malloc(N*sizeof(float));
    
    
    // this is for holding onto overlap portion each iteration
    float* overlapBuffer;
    overlapBuffer = (float*)malloc((N/2)*sizeof(float));
    
    // generic scratch buffer
    float* scratchBuffer;
    scratchBuffer = (float*)malloc(N*sizeof(float));
    
    // these 1024 samples get compacted into 512 reals and imaginaries
    DSPSplitComplex holdingBufferSplit;
    holdingBufferSplit.realp = new float[N/2];
    holdingBufferSplit.imagp = new float[N/2];
    
    // FFT-ed filter for overlap save (need N = 1024 length)
    // https://stackoverflow.com/questions/2929401/dsp-filtering-in-the-frequency-domain-via-fft
        
    // writeOutput(filterCoeffs, PADDED_LENGTH, 44100);
    
    DSPSplitComplex filterSplit;
    filterSplit.realp = new float[N/2];
    filterSplit.imagp = new float[N/2];
    vDSP_ctoz((DSPComplex *)filterCoeffs, 2, &filterSplit, 1, N/2);
    vDSP_fft_zrip(setup, &filterSplit, 1, log2n, FFT_FORWARD);
    
    
    
    // sanity check that filter magnitude makes sense
    /*
    float* filterMag;
    filterMag = new float[N];
    for(int i=0; i<N/2; i++){
        filterMag[i] = fabsf(filterSplit.realp[i]*filterSplit.realp[i] + filterSplit.imagp[i]*filterSplit.imagp[i]);
    }
    // N/2 samples is the positive frequency section of the filter magnitude response plot
    writeOutput(filterMag, N/2, 44100);
    */
    
    // 1024 length
    float* reInterleaved;
    reInterleaved = new float[N];
    
    [self.audioManager setOutputBlock:^(float *outData,
                                        UInt32 numFrames,
                                        UInt32 numChannels) {
        
        

        
        
        
        
        // 2. process the input ring buffer, producing samples for the OUTPUT ring buffer
        
        // while the input ring buffer contains enough to run an N==1024 length FFT
        while(wself.ringBufferIn->NumUnreadFrames() >= N) {
                        
            // a. run the FFT on the (windowed!) samples at the end of the input ring buffer
            
            // fetch head is advancing 512 samples each time so it's continguous with itself
            wself.ringBufferIn->FetchData(bufferInL, numFrames, 0, 1);
            
            
            for(int i=0; i<N; i++) {
                if(i<N/2);
//                    bufferInL[i] *= myWindow[i];
                else
                    bufferInL[i] = 0.0;
            }
            
            // sanity check that time domain output is correct
            
//            if(!dumpHit){
//                if(dumpCount==dumpMatch) {
//                    dumpHit = 1;
//                    writeOutput(bufferInL, 2*numFrames, 44100);
//                }
//                else
//                    dumpCount++;
//            }
            
            
            // pack this into format expected by upcoming FFT
            vDSP_ctoz((DSPComplex*)bufferInL, 2,
                      &holdingBufferSplit, 1,
                      N/2);

            
            // do the FFT. data is now in holdingBufferSplit (real & imag)
            vDSP_fft_zrip(setup,
                          &holdingBufferSplit, 1,
                          log2n,
                          FFT_FORWARD);
            
            // sanity check that frequency domain signal is correct

//            if(!dumpHit){
//                if(dumpCount==dumpMatch) {
//                    dumpHit = 1;
//                    for(int i=0; i<N/2; i++){
//                        scratchBuffer[i] = fabsf(
//                                                 holdingBufferSplit.realp[i]*holdingBufferSplit.realp[i] +
//                                                 holdingBufferSplit.imagp[i]*holdingBufferSplit.imagp[i]
//                                                 );
//                    }
//                    writeOutput(scratchBuffer, N, 44100);
//                }
//                else
//                    dumpCount++;
//            }
            

            
            float preserveSigNyq = holdingBufferSplit.imagp[0];
            holdingBufferSplit.imagp[0] = 0.0f;

            float preserveFiltNyq = filterSplit.imagp[0];
            filterSplit.imagp[0] = 0.0f;
            
            // do FFT(x) * FFT(h) with vmul
//            vDSP_zvmul(&holdingBufferSplit, 1,
//                       &filterSplit, 1,
//                       &holdingBufferSplit, 1,
//                       N/2,
//                       1);
            
            holdingBufferSplit.imagp[0] = preserveFiltNyq * preserveSigNyq;
            filterSplit.imagp[0] = preserveFiltNyq;
            
            // sanity check that frequency domain product of signal and filter is correct
            
//            if(!dumpHit){
//                if(dumpCount==dumpMatch) {
//                    dumpHit = 1;
//                    for(int i=0; i<N/2; i++){
//                        scratchBuffer[i] = fabsf(
//                                                      holdingBufferSplit.realp[i]*holdingBufferSplit.realp[i] +
//                                                      holdingBufferSplit.imagp[i]*holdingBufferSplit.imagp[i]
//                                                      );
//                    }
//                    writeOutput(scratchBuffer, N, 44100);
//                }
//                else
//                    dumpCount++;
//            }
            
            // inverse FFT the product
            vDSP_fft_zrip(setup,
                          &holdingBufferSplit, 1,
                          log2n,
                          FFT_INVERSE);
            
            // unpack the IFFT data
            vDSP_ztoc(&holdingBufferSplit, 1,
                      (DSPComplex*)bufferOutL, 2,
                      N/2);
            
            // may need to play around with this scaling factor
            float scale = 0.25/N;
            vDSP_vsmul(bufferOutL, 1,
                       &scale,
                       bufferOutL, 1,
                       N);
            
            // sanity check that inversed time domain samples are correct
            
//            if(!dumpHit){
//                if(dumpCount==dumpMatch) {
//                    dumpHit = 1;
//                    writeOutput(bufferOutL, N, 44100);
//                }
//                else
//                    dumpCount++;
//            }
            
            
            
            
            
            // at this point, each IFFT is FILTER_LENGTH/2 delayed
            // de-offset the IFFT'd time domain samples by FILTER_LENGTH/2
//            for (int i=0; i<N; i++){
//                bufferOutL[i] = bufferOutL[ (i + FILTER_LENGTH/2) % N ];
//            }
            
            
//            if(!dumpHit){
//                if(dumpCount==dumpMatch) {
//                    dumpHit = 1;
//                    writeOutput(bufferOutL, N, 44100);
//                }
//                else
//                    dumpCount++;
//            }
            
            
            
            
            
            
            
            
            // do the overlap add
            // grab 2nd half from last iteration's output and put it into the overlap buffer
            wself.ringBufferOverlap->FetchData(overlapBuffer, numFrames, 0, 1);
            
//            if(!dumpHit){
//                if(dumpCount==dumpMatch) {
//                    dumpHit = 1;
//                    writeOutput(overlapBuffer, numFrames, 44100);
//                }
//                else
//                    dumpCount++;
//            }
            
            
            // current 1st half output PLUS old 2nd half output
            for(int i=0; i<numFrames; i++){
                bufferOutL[i] += overlapBuffer[i];
            }
            
            // put the overlap-add data into the output ring buffer
            // have to reinterleave it though
            for(int i=0; i<numFrames; i++){
                reInterleaved[2*i] = bufferOutL[i];
                reInterleaved[2*i+1] = bufferOutL[i];
            }
            wself.ringBufferOut->AddNewInterleavedFloatData(reInterleaved, numFrames, numChannels);
            
//            if(!dumpHit){
//                if(dumpCount==dumpMatch) {
//                    dumpHit = 1;
//                    writeOutput(bufferOutL, N/2, 44100);
//                }
//                else
//                    dumpCount++;
//            }
            
            // need to store the 2nd half of current bufferOutL for next time
            for(int i=0; i<numFrames; i++){
                overlapBuffer[i] = bufferOutL[i+N/2];
            }
            wself.ringBufferOverlap->AddNewFloatData(overlapBuffer, numFrames, 0);
            
//            if(!dumpHit){
//                if(dumpCount==dumpMatch) {
//                    dumpHit = 1;
//                    writeOutput(overlapBuffer, N/2, 44100);
//                }
//                else
//                    dumpCount++;
//            }
            
        }
        
        // this block is for testing.
        // comment all the output related stuff above before uncommenting and running this
        /*
        // reinterleave it
        wself.ringBufferIn->FetchData(deInterleavedL, numFrames, 0, 1);
        wself.ringBufferIn->FetchData(deInterleavedR, numFrames, 1, 1);
        
        
        
        for(int i=0; i<numFrames; i++){
            reInterleaved[2*i] = deInterleavedL[i];
            reInterleaved[2*i+1] = deInterleavedR[i];
        }
        
        if(!dumpHit){
            if(dumpCount==dumpMatch) {
                dumpHit = 1;
                writeOutput(reInterleaved, N, 44100);
            }
            else
                dumpCount++;
        }
        
        wself.ringBufferOut->AddNewInterleavedFloatData(reInterleaved, numFrames, numChannels);
         */
        
        
        
        wself.ringBufferOut->FetchInterleavedData(outData, numFrames, numChannels);
         
        
        
        
        
        
        
    }];
    
    
 

    // START IT UP YO
    [self.audioManager play];
}

@end
