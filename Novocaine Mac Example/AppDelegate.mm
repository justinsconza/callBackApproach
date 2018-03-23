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
    if (_ringBuffer){
        delete _ringBuffer;
    }
}

- (void)applicationDidFinishLaunching:(NSNotification *)aNotification
{
    
    self.audioManager = [Novocaine audioManager];
    
//    self.ringBuffer = new RingBuffer(64*512, 2);
    self.ringBuffer = new RingBuffer(16*512, 2);
    
    __weak AppDelegate * wself = self;
    
    typedef NS_ENUM(NSUInteger,DspAlgorithm) {
        basic,
        delay,
        fftIfft
    };
    
    
    DspAlgorithm myDspAlgorithm = fftIfft;
    
    // --------------------------------------------------------------- //
    // ------------------ BASIC THROUGHPUT --------------------------- //
    // --------------------------------------------------------------- //
    if (myDspAlgorithm == basic) {
        // Basic playthru example
        [self.audioManager setInputBlock:^(float *data,
                                           UInt32 numFrames,
                                           UInt32 numChannels) {
            float volume = 0.5;
            vDSP_vsmul(data, 1, &volume, data, 1, numFrames*numChannels);
            wself.ringBuffer->AddNewInterleavedFloatData(data, numFrames, numChannels);
        }];

        [self.audioManager setOutputBlock:^(float *outData,
                                            UInt32 numFrames,
                                            UInt32 numChannels) {
            wself.ringBuffer->FetchInterleavedData(outData, numFrames, numChannels);
        }];
    }
    
    // --------------------------------------------------------------- //
    // ------------------ SLAPBACK DELAY ----------------------------- //
    // --------------------------------------------------------------- //
    
    if (myDspAlgorithm == delay) {
        
        // A simple delay that's hard to express without ring buffers
        // numFrames is 512 and numChannels is 2
        // data points to the input buffer used by the input callback
        [self.audioManager setInputBlock:^(float *data, UInt32 numFrames, UInt32 numChannels) {
            
            // AddNewInterleavedFloatData moves data from input buffer called data to the ringBuffer
            wself.ringBuffer->AddNewInterleavedFloatData(data, numFrames, numChannels);
            
        }];
    
        int echoDelay = 11025;
        float *holdingBuffer = (float *)calloc(16384, sizeof(float));
    
        // outData points to the output buffer used by the render callback
        [self.audioManager setOutputBlock:^(float *outData, UInt32 numFrames, UInt32 numChannels) {
            
            // would like to know how separated these heads actually are
            // put a printf in call to add new data and print mNumUnreadFrames[0] with %lld
            
            /*
            step 1: grab a block of 512 samples starting from current idx in ring buffer.
                    put this block into the output buffer.
            */
            
            // this wasn't calling with numChannels, so had to add that.  CRITICAL!
            // FetchInterleaveData moves data from ringBuffer to output buffer called outData
            wself.ringBuffer->FetchInterleavedData(outData, numFrames, numChannels);
            
            float volume = 0.8;
            vDSP_vsmul(outData, 1, &volume, outData, 1, numFrames*numChannels);
        
            /*
            step 2: go back from last written index by one delay length
                    and one block length (512) and place the read head there
            */
            for (int i=0; i<numChannels; i++){
                wself.ringBuffer->SeekReadHeadPosition(-echoDelay-numFrames, i);
                
            }
            
            /*
            step 3: grab a block of 512 samples starting where read head was just placed.
                    put this block into the holding buffer.
            */
            wself.ringBuffer->FetchInterleavedData(holdingBuffer, numFrames, numChannels);
            
            /*
            step 4: the last written index is now one delay length behind
                    the idx we started step 1 with. so move the read head
                    to one position ahead of where we started out in the beginning of step 1
            */
            for (int i=0; i<numChannels; i++){
                wself.ringBuffer->SeekReadHeadPosition(echoDelay,i);
            }
            
            // turn down the holding buffer samples by 1/2
            volume = 0.5;
            vDSP_vsmul(holdingBuffer, 1, &volume, holdingBuffer, 1, numFrames*numChannels);
            
            // add the holding buffer block to the output block.  this will get picked up by the renderCallback and played out by the speakers.
            vDSP_vadd(holdingBuffer, 1, outData, 1, outData, 1, numFrames*numChannels);
            
        }];
    }
    
    // --------------------------------------------------------------- //
    // ------------------ FFT TO IFFT -------------------------------- //
    // --------------------------------------------------------------- //
    
    if (myDspAlgorithm == fftIfft) {
        
        [self.audioManager setInputBlock:^(float *data,
                                           UInt32 numFrames,
                                           UInt32 numChannels) {
            float volume = 0.5;
            vDSP_vsmul(data, 1, &volume, data, 1, numFrames*numChannels);
            wself.ringBuffer->AddNewInterleavedFloatData(data, numFrames, numChannels);
        }];
        
        int N = 1024;
        int log2n = 10;    // 2^9 = 512
        
        DSPSplitComplex zSplit;
        zSplit.realp = new float[N/2];
        zSplit.imagp = new float[N/2];
        
        float* mag;
        float* phase;
        mag = new float[N/2];
        phase = new float[N/2];
        
    
        float* myWindow;
        myWindow = new float[N];
        vDSP_hann_window(myWindow, N, vDSP_HANN_DENORM);
        

        // prepare the fft algo (you want to reuse the setup across fft calculations)
        FFTSetup setup = vDSP_create_fftsetup(log2n, kFFTRadix2);
    
        [self.audioManager setOutputBlock:^(float *outData,
                                            UInt32 numFrames,
                                            UInt32 numChannels) {
            
            
            // copy 512 time domain samples to output buffer
            wself.ringBuffer->FetchInterleavedData(outData, numFrames, numChannels);
            
            for (int iChannel = 0; iChannel<numChannels; iChannel++){
                
//                vDSP_vmul(outData, 1, myWindow, 1, outData + iChannel, 1, N);
//                for (int i=0; i<N; i++){
//                    (outData + iChannel)[i] *= myWindow[i];
//                }
                
                // get it into packed format for fft
                vDSP_ctoz((DSPComplex *) (outData + iChannel), 2, &zSplit, 1, N/2);
                
                // calculate the fft
                vDSP_fft_zrip(setup, &zSplit, 1, log2n, FFT_FORWARD);
                
                
                
                for(int i=0; i<N/2; i++){
                    
                    mag[i] = sqrtf(zSplit.realp[i]*zSplit.realp[i] + zSplit.imagp[i]*zSplit.imagp[i]);
                    phase[i] = atan2f(zSplit.imagp[i], zSplit.realp[i]);
                     
                }
                
                
                float* zReal = zSplit.realp;
                float* zImag = zSplit.imagp;
                
                for(int i=0; i<N/2; i++){
                    *zReal++ = mag[i]*cosf(phase[i]);
                    *zImag++ = mag[i]*sinf(phase[i]);
                }
                
                
                // calculate the ifft. leaves it in packed format.
                vDSP_fft_zrip(setup, &zSplit, 1, log2n, FFT_INVERSE);
                
                // unpack it into a real vector
                vDSP_ztoc(&zSplit, 1, (DSPComplex *) (outData + iChannel), 2, N/2);
                
                // scale the result
                float scale = 0.5/N;
                vDSP_vsmul(outData + iChannel, 1, &scale, outData + iChannel, 1, N);

            }
            
        }];
        
    }
 
    // AUDIO FILE READING COOL!
    // ========================================    
//    NSURL *inputFileURL = [[NSBundle mainBundle] URLForResource:@"TLC" withExtension:@"mp3"];
//
//    self.fileReader = [[AudioFileReader alloc]
//                       initWithAudioFileURL:inputFileURL
//                       samplingRate:self.audioManager.samplingRate
//                       numChannels:self.audioManager.numOutputChannels];
//
//    self.fileReader.currentTime = 5;
//    [self.fileReader play];
//
//
//    __block int counter = 0;
//
//
//    [self.audioManager setOutputBlock:^(float *data, UInt32 numFrames, UInt32 numChannels)
//     {
//         [wself.fileReader retrieveFreshAudio:data numFrames:numFrames numChannels:numChannels];
//         counter++;
//         if (counter % 80 == 0)
//             NSLog(@"Time: %f", wself.fileReader.currentTime);
//
//     }];
    
    
    // AUDIO FILE WRITING YEAH!
    // ========================================    
//    NSArray *pathComponents = [NSArray arrayWithObjects:
//                               [NSSearchPathForDirectoriesInDomains(NSDocumentDirectory, NSUserDomainMask, YES) lastObject], 
//                               @"My Recording.m4a", 
//                               nil];
//    NSURL *outputFileURL = [NSURL fileURLWithPathComponents:pathComponents];
//
//    self.fileWriter = [[AudioFileWriter alloc]
//                       initWithAudioFileURL:outputFileURL 
//                       samplingRate:self.audioManager.samplingRate
//                       numChannels:self.audioManager.numInputChannels];
//    
//    
//    __block int counter = 0;
//    self.audioManager.inputBlock = ^(float *data, UInt32 numFrames, UInt32 numChannels) {
//        [wself.fileWriter writeNewAudio:data numFrames:numFrames numChannels:numChannels];
//        counter += 1;
//        if (counter > 10 * wself.audioManager.samplingRate / numChannels) { // 10 seconds of recording
//            wself.audioManager.inputBlock = nil;
//        }
//    };

    // START IT UP YO
    [self.audioManager play];
}

@end
