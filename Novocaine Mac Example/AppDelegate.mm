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
    
    int ringBufferLength = 32*512;
    self.ringBufferIn = new RingBuffer(ringBufferLength, 2);
    self.ringBufferOut = new RingBuffer(ringBufferLength, 2);
    self.ringBufferOverlap = new RingBuffer(512, 2);
    
    __weak AppDelegate * wself = self;
 
   
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
    for(int i=0; i<N; i++){
        deInterleavedL[i] = deInterleavedR[i] = 0.0;
    }
    
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
        
        // wself.ringBufferIn->AddNewInterleavedFloatData(data, numFrames, numChannels);
        
    }];
  
    


    // don't think i need a window
    float* myWindow;
    myWindow = new float[N/2];
    vDSP_hann_window(myWindow, N/2, vDSP_HANN_DENORM);
    
    FFTSetup setup = vDSP_create_fftsetup(log2n, kFFTRadix2);
    
    // doing overlap add, grab 512 new and concatenate with 512 old
    float* bufferInL;
    bufferInL = (float*)malloc(N*sizeof(float));
    float* bufferInR;
    bufferInL = (float*)malloc(N*sizeof(float));
    
    float* bufferOutL;
    bufferOutL = (float*)malloc(N*sizeof(float));
    float* bufferOutR;
    bufferOutR = (float*)malloc(N*sizeof(float));
    
    // this is for holding onto overlap portion each iteration
    float* overlapBuffer;
    overlapBuffer = (float*)malloc((N/2)*sizeof(float));
    
    // generic 1024 long scratch buffer
    float* scratchBuffer;
    scratchBuffer = (float*)malloc(N*sizeof(float));

    // 1024 length
    float* reInterleaved;
    reInterleaved = new float[N];
    
    for(int i=0; i<N; i++){
        bufferInL[i] = bufferInL[i] = 0.0;
        bufferOutL[i] = bufferOutR[i] = 0.0;
        scratchBuffer[i] = 0.0;
        reInterleaved[i] = 0.0;
        if (i<N/2)
            overlapBuffer[i] = 0.0;
    }
    
    

   
    // these 1024 samples get compacted into 512 reals and imaginaries
    DSPSplitComplex holdingBufferSplit;
    holdingBufferSplit.realp = new float[N/2];
    holdingBufferSplit.imagp = new float[N/2];
        
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
    

    __block int counter = 10;
    __block int doneCounting = 0;
    
    [self.audioManager setOutputBlock:^(float *outData,
                                        UInt32 numFrames,
                                        UInt32 numChannels) {
        
        
        // need to avoid a pop sound at the beginning
        if(!doneCounting) {
            counter -= 1;
            if(counter == 0)
                doneCounting = 1;
        }
        
        if(doneCounting) {
        
        
        
            // NEED 512 NEW SAMPLES TO RENDER AN OUTPUT BLOCK
            while(wself.ringBufferIn->NumUnreadFrames() >= N/2) {
                
                // ---------- GET INPUT CORRECT ---------------------------- //
                
                // FETCH 512 SAMPLES FROM INPUT RING BUFFER --> bufferInL
                // printf("take from input: %d\n", wself.ringBufferIn->ReadHeadPosition(0));
                wself.ringBufferIn->FetchData(bufferInL, numFrames, 0, 1);
                

                // FETCH OLD 512 SAMPLES FROM LAST TIME --> overlapBuffer
                // printf("begin reading at: %d\n", wself.ringBufferOverlap->ReadHeadPosition(0));
                wself.ringBufferOverlap->FetchData(overlapBuffer, numFrames, 0, 1);
                // printf("finished reading at: %d\n", wself.ringBufferOverlap->ReadHeadPosition(0));
                
                // 2ND HALF IS NEW 512, 1ST HALF IS OLD 512
                for(int i=0; i<numFrames; i++) {
                    bufferInL[i+numFrames] = bufferInL[i];
                    bufferInL[i] = overlapBuffer[i];
                }
                
    //            if(dumpCount++ == dumpMatch)
    //                writeOutput(bufferInL, N, 44100);
                
                // --------- FREQUENCY DOMAIN SETUP ------------------------ //
                
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
    //            if(dumpCount++ == dumpMatch) {
    //                for(int i=0; i<N/2; i++){
    //                    scratchBuffer[i] = fabsf(
    //                                             holdingBufferSplit.realp[i]*holdingBufferSplit.realp[i] +
    //                                             holdingBufferSplit.imagp[i]*holdingBufferSplit.imagp[i]
    //                                             );
    //                }
    //                writeOutput(scratchBuffer, N, 44100);
    //            }
                
                // --------- FREQUENCY DOMAIN MULTIPLY --------------------- //
                
                float preserveSigNyq = holdingBufferSplit.imagp[0];
                holdingBufferSplit.imagp[0] = 0.0f;
                
                float preserveFiltNyq = filterSplit.imagp[0];
                filterSplit.imagp[0] = 0.0f;
                
     
                // do FFT(x) * FFT(h) with vmul
                vDSP_zvmul(&holdingBufferSplit, 1,
                           &filterSplit, 1,
                           &holdingBufferSplit, 1,
                           N/2,
                           1);
                
                holdingBufferSplit.imagp[0] = preserveFiltNyq * preserveSigNyq;
                filterSplit.imagp[0] = preserveFiltNyq;
                
                
                // ---------------------------------------------------- //
                
                // sanity check that frequency domain product of signal and filter is correct
    //            if(dumpCount++ == dumpMatch) {
    //                for(int i=0; i<N/2; i++){
    //                    scratchBuffer[i] = fabsf(
    //                                             holdingBufferSplit.realp[i]*holdingBufferSplit.realp[i] +
    //                                             holdingBufferSplit.imagp[i]*holdingBufferSplit.imagp[i]
    //                                             );
    //                }
    //                writeOutput(scratchBuffer, N, 44100);
    //            }

                
                // inverse FFT the product
                vDSP_fft_zrip(setup,
                              &holdingBufferSplit, 1,
                              log2n,
                              FFT_INVERSE);
                
                // unpack the IFFT data
                vDSP_ztoc(&holdingBufferSplit, 1,
                          (DSPComplex*)scratchBuffer, 2,
                          N/2);
                
                // may need to play around with this scaling factor
                float scale = 0.25/N;
                vDSP_vsmul(scratchBuffer, 1,
                           &scale,
                           scratchBuffer, 1,
                           N);                  // should this be N/2 ?????

                // --------- PREPARE FOR NEXT TIME & THE OUTPUT ------------ //
                
                // need to store the 2nd half of current input for next time
                // also save 2nd half of output (doing overlap-save) / reinterleave
                for(int i=0; i<numFrames; i++){
                    overlapBuffer[i] = bufferInL[i+N/2];
                    reInterleaved[2*i] = scratchBuffer[i+N/2];
                    reInterleaved[2*i+1] = scratchBuffer[i+N/2];
    //                bufferOutL[i] = scratchBuffer[i+N/2];
                }
                
    //            if(dumpCount++ == dumpMatch)
    //                writeOutput(bufferOutL, N/2, 44100);
                
                // printf("begin writing at: %d\n", wself.ringBufferOverlap->WriteHeadPosition(0));
                wself.ringBufferOverlap->AddNewFloatData(overlapBuffer, numFrames, 0);
                // printf("end writing at: %d\n", wself.ringBufferOverlap->WriteHeadPosition(0));

                wself.ringBufferOut->AddNewInterleavedFloatData(reInterleaved, N/2, numChannels);

            } // end while 512 needed sample while loop
            
        } // end de-pop conditional
        
//         float volume = 1.0;
//         vDSP_vsmul(outData, 1, &volume, outData, 1, numFrames*numChannels);
     
        wself.ringBufferOut->FetchInterleavedData(outData, numFrames, numChannels);
      
    }];

    // START IT UP YO
    [self.audioManager play];
}

@end

