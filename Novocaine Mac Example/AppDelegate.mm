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
    self.ringBufferBetween = new RingBuffer(4*512, 2);
    
    __weak AppDelegate * wself = self;
 
   

    
    
    // numFrames is 512, numChannels is 2
        
    [self.audioManager setInputBlock:^(float *data,
                                       UInt32 numFrames,
                                       UInt32 numChannels) {
        
        // 1. stash all input samples into an INPUT ring buffer
        /*
        float volume = 0.75;
        vDSP_vsmul(data, 1, &volume, data, 1, numFrames*numChannels);
        */
        wself.ringBufferIn->AddNewInterleavedFloatData(data, numFrames, numChannels);
        
    }];
  
    
    // doing FFT's that are twice as long as numFrames == 512
    
    int N = 1024;       // 2N is 2*numFrames == 2*512
    int log2n = 10;     // 2^10 = 1024

    // don't think i need a window but if i do, it should be same length as FFT's == 1024
    float* myWindow;
    myWindow = new float[N];
    vDSP_hann_window(myWindow, N, vDSP_HANN_DENORM);
    
    FFTSetup setup = vDSP_create_fftsetup(log2n, kFFTRadix2);
    
    // doing overlap save, need to grab 1024 time samples each go round
    float* holdingBufferIn;
    holdingBufferIn = (float*)malloc(N*sizeof(float));
    
    float* holdingBufferOut;
    holdingBufferOut = (float*)malloc(N*sizeof(float));
    
    // zero out second half of holdingBuffer for fake linear convolution
    /*
    for(int i=N/2; i<N; i++){
        holdingBuffer[i] = 0.0;
    }
    */
    
    // this is for holding onto 512 input time samples for overlap add remembering each iteration
    float* scratchBuffer;
    scratchBuffer = (float*)malloc((N/2)*sizeof(float));
    
    // this is for verifying overlap add stuff is continguous across iterations
    float* scratchBuffer2;
    scratchBuffer2 = (float*)malloc(N*sizeof(float));
    
    // these 1024 samples get compacted into 512 reals and imaginaries
    DSPSplitComplex holdingBufferSplit;
    holdingBufferSplit.realp = new float[N/2];
    holdingBufferSplit.imagp = new float[N/2];
    
    // FFT-ed filter for overlap save (need N = 2*512 length)
    // https://stackoverflow.com/questions/2929401/dsp-filtering-in-the-frequency-domain-via-fft
        
    // writeOutput(filterCoeffs, N, 44100);
    
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
    
    __block int dumpHit = 0;
    __block int dumpMatch = 10;
    __block int dumpCount = 1;
    
    [self.audioManager setOutputBlock:^(float *outData,
                                        UInt32 numFrames,
                                        UInt32 numChannels) {
        
        
        
        // 2. process the input ring buffer, producing samples for the OUTPUT ring buffer
        
        // while the input ring buffer contains enough to run an N==1024 length FFT
        while(wself.ringBufferIn->NumUnreadFrames() >= N){
                        
            // a. run the FFT on the (windowed!) samples at the end of the input ring buffer
            
            // 1024 new input samples, first 512 are overlapped with previous iterations last 512
            // rewind the read head N/2 because reading N and have to overlap each read with previous iteration
            wself.ringBufferIn->SeekReadHeadPosition(-numFrames);
            wself.ringBufferIn->FetchData(holdingBufferIn, 2*numFrames, 0, 1);
            
            // don't know if i need to window the signal or not...
            // DO NOT WINDOW
            /*
            vDSP_vmul(holdingBufferIn, 1,
                      myWindow, 1,
                      holdingBufferIn, 1,
                      N);
            */
            
            // sanity check that time domain output is correct
            
//            if(!dumpHit){
//                if(dumpCount==dumpMatch) {
//                    dumpHit = 1;
//                    writeOutput(holdingBufferIn, N, 44100);
//                }
//                else
//                    dumpCount++;
//            }
            
            
            // pack this into format expected by upcoming FFT
            vDSP_ctoz((DSPComplex*)holdingBufferIn, 2,
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
            
            
            
            // do FFT(x) * FFT(h) with vmul
//            vDSP_zvmul(&holdingBufferSplit, 1,
//                       &filterSplit, 1,
//                       &holdingBufferSplit, 1,
//                       N/2,
//                       1);
            
            
            // sanity check that frequency domain product of signal and filter is correct
            
//            if(!dumpHit){
//                if(dumpCount==dumpMatch) {
//                    dumpHit = 1;
//                    for(int i=0; i<N/2; i++){
//                        printf("%f\n",holdingBufferSplit.realp[i]);
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
                      (DSPComplex*)holdingBufferOut, 2,
                      N/2);
            
            // may need to play around with this scaling factor
            float scale = 0.25/N;
            vDSP_vsmul(holdingBufferOut, 1,
                       &scale,
                       holdingBufferOut, 1,
                       N);
            
            // sanity check that inversed time domain samples are correct
            
//            if(!dumpHit){
//                if(dumpCount==dumpMatch) {
//                    dumpHit = 1;
//                    writeOutput(holdingBufferOut, N, 44100);
//                }
//                else
//                    dumpCount++;
//            }
            
            
            // at this point, each IFFT is FILTER_LENGTH/2 delayed
            
            
            // de-offset the IFFT'd time domain samples by FILTER_LENGTH/2
            // OFFSET IS WRONG DON'T DO IT
            /*
            for (int i=0; i<N; i++){
                holdingBufferOut[i] = holdingBufferOut[ (i + FILTER_LENGTH/2) % N ];
            }
            */
            
//            if(!dumpHit){
//                if(dumpCount==dumpMatch) {
//                    dumpHit = 1;
//                    writeOutput(holdingBufferOut, N, 44100);
//                }
//                else
//                    dumpCount++;
//            }
            
            
            // do the overlap add
            // grab 2nd half from last iteration and put it into the scratch buffer
            // LEARNED OFFSET IS BAD FROM THIS CODE HERE
            wself.ringBufferBetween->FetchData(scratchBuffer, numFrames, 0, 1);
            
//            if(!dumpHit){
//                if(dumpCount==dumpMatch) {
//                    dumpHit = 1;
//                    writeOutput(scratchBuffer, N/2, 44100);
//                }
//                else
//                    dumpCount++;
//            }
            
            
            // current 1st half (holdingBufferOut 1st half) + previous 2nd half (scratchBuffer)
            for(int i=0; i<N/2; i++){
                scratchBuffer[i] += holdingBufferOut[i];
            }
            
            // put the overlap-add data from the scratch buffer into the output ring buffer
            wself.ringBufferOut->AddNewFloatData(scratchBuffer, numFrames, 0);
            
//            if(!dumpHit){
//                if(dumpCount==dumpMatch) {
//                    dumpHit = 1;
//                    writeOutput(scratchBuffer, N/2, 44100);
//                }
//                else
//                    dumpCount++;
//            }
            
            // need to store the 2nd half of current holdingBufferOut for next time
            for(int i=0; i<N/2; i++){
                scratchBuffer[i] = holdingBufferOut[i+N/2];
            }
            wself.ringBufferBetween->AddNewFloatData(scratchBuffer, numFrames, 0);
            
//            if(!dumpHit){
//                if(dumpCount==dumpMatch) {
//                    dumpHit = 1;
//                    writeOutput(scratchBuffer, N/2, 44100);
//                }
//                else
//                    dumpCount++;
//            }
            
        }

        // finally, put the output ring buffer data into the output stream
        wself.ringBufferOut->FetchInterleavedData(outData, numFrames, numChannels);
        
//        if(!dumpHit){
//            if(dumpCount==dumpMatch) {
//                dumpHit = 1;
//                writeOutput(outData, N/2, 44100);
//            }
//            else
//                dumpCount++;
//        }
        
        
    }];
    
    
 

    // START IT UP YO
    [self.audioManager play];
}

@end
