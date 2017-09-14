#include <complex.h>
#include <fftw3.h>
#include <stdlib.h>
#include <stdio.h>
#include <sndfile.h>
#include <portaudio.h>

#include "chebyshev.h"

// Struct to hold data for audio callback
typedef struct {
    float *data;
    size_t len;
    size_t ptr;
} Data;

// Audio callback
int audio_callback(const void *input,
                   void *output,
                   unsigned long frame_count,
                   const PaStreamCallbackTimeInfo *timeInfo,
                   PaStreamCallbackFlags statusFlags,
                   void *userData)
{
    /* Cast data passed through stream to our structure. */
    Data *data = (Data *)userData;
    float *out = (float *)output;
    unsigned int i;
    (void) input; /* Prevent unused variable warning */

    for (i = 0; i < frame_count; ++i)
    {
        if (data->ptr < data->len)
            out[i] = data->data[data->ptr++];
        else
            out[i] = 0;
    }
    if (data->ptr < data->len)
        return 0;
    else
        return paComplete;
}

// Normalize an fft (only necessary when going forwards, then backwards)
void normalize(fftw_complex *a, int n)
{
    for (int i = 0; i < n; i++)
    {
        a[i] /= n;
    }
}

// Prints the port audio error corresponding to input parameter
void print_pa_error(PaError err)
{
    fprintf(stderr, "PortAudio error: %s\n", Pa_GetErrorText(err));
}

int main(void)
{
    
    // Load in test wav file
    SF_INFO sfinfo = {0};
    SNDFILE *f = sf_open("test.wav", SFM_READ, &sfinfo);

    // Check if it loaded
    if (f == NULL)
    {
        fprintf(stderr, "ERROR: %s\n\n", sf_strerror(f));
        return -1;
    }

    // Retrieve samples from file
    int samples = sfinfo.frames*sfinfo.channels;
    float frames[samples];
    sf_count_t len_read = sf_readf_float(f, frames, sfinfo.frames);

    // Check if whole file was read
    if (len_read < sfinfo.frames)
    {
        fprintf(stderr, "ERROR: Could not read all frames!\n");
    }

    // Length of audio in samples
    int n = sfinfo.frames;

    // Variables used for FFT
    fftw_complex *in, *out;
    fftw_plan pforward;

    // FFT length
    int fft_len = 1024;

    // Create FFT arrays and FFT plan
    in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*fft_len);
    out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*fft_len);
    pforward = fftw_plan_dft_1d(fft_len, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Sum stereo channels down to mono
    double *song = malloc(sizeof(double)*n);
    for (int i = 0; i < n; ++i)
        song[i] = (frames[2*i] + frames[2*i+1])/2;

    // Set filter parameters
    int sample_rate = sfinfo.samplerate; 
    double cutoff = 1000;                   // Cutoff frequency 
    double pr = 0.5;                        // Percent ripple (0.5%)
    int np = 6;                             // Number of poles
    // Get chebyshev coefficients
    FLT_COEFF ret = chebyshev((float) cutoff / sample_rate, false, pr, np);

    // Create array to hold filtered song
    const int filtered_len = n+np;
    double *filtered = malloc(sizeof(double)*filtered_len);
    for (int i = 0; i < filtered_len; ++i)
        filtered[i] = 0;

    // Apply filter
    apply_filter(song, n, ret.a, ret.b, np+1, filtered, filtered_len);

    // Put samples in FFT input array
    for (int i = 0 ; i < fft_len; ++i)
    {
        in[i] = song[i];
    }

    // Run FFT for unfiltered song
    fftw_execute(pforward);
    //normalize(out, fft_len);
    
    // Open output data file
    FILE *out_f = fopen("data.dat", "w");
    
    // Print sample rate to convert fft bin numbers to frequencies
    fprintf(out_f, "%d\n", sample_rate);

    // Output magnitude of each fft bin
    for (int i = 0; i < fft_len; ++i)
        fprintf(out_f, "%f ", cabs(out[i]));

    // Next line
    fprintf(out_f, "\n");

    // Put filtered samples in FFT input array
    for (int i = 0; i < fft_len; ++i)
    {
        in[i] = filtered[i];
    }

    // Run FFT for filtered song
    fftw_execute(pforward);
    //normalize(out, fft_len);

    // Output magnitude of each fft bin
    for (int i = 0; i < fft_len; ++i)
        fprintf(out_f, "%f ", cabs(out[i]));

    // Close output file
    fclose(out_f);

    // Create buffer to hold song, some silence, and the filtered song
    Data data;
    const size_t num_zero = sfinfo.samplerate;
    data.len = n + filtered_len + num_zero;
    data.ptr = 0;
    data.data = malloc(sizeof(float)*(n + filtered_len + num_zero));
    for (int i = 0; i < n; ++i)
        data.data[i] = song[i];
    for (int i = 0; i < num_zero; ++i)
        data.data[n+i] = 0;
    for (int i = 0; i < filtered_len; ++i)
        data.data[n+num_zero+i] = filtered[i];


    // Initialize Port Audio to play audio
    PaError err = Pa_Initialize();
    if( err != paNoError){
        print_pa_error(err);
        goto error;
    }

    // Open audio stream
    PaStream *stream;
    err = Pa_OpenDefaultStream(&stream,
                               0,                               // # input channels
                               1,                               // # output channels
                               paFloat32,                       // Type
                               sfinfo.samplerate,
                               paFramesPerBufferUnspecified,    // Frames per buffer
                               audio_callback,                  // Callback function for audio
                               &data);                          // data to be passed to callback
    if( err != paNoError){
        print_pa_error(err);
        goto error;
    }

    // Start audio
    err = Pa_StartStream(stream);
    if( err != paNoError){
        print_pa_error(err);
        goto error;
    }

    // Wait until audio is done playing
    while (Pa_IsStreamActive(stream) == 1);

    // Close stream
    err = Pa_CloseStream(stream);
    if( err != paNoError){
        print_pa_error(err);
        goto error;
    }

    // Terminate Port Audio session
    err = Pa_Terminate();
    if( err != paNoError){
        print_pa_error(err);
        goto error;
    }

error:

    // Free all resources
    free(song);
    free(filtered);
    free(ret.a);
    free(ret.b);
    free(data.data);

    sf_close(f);
    fftw_destroy_plan(pforward); 
    fftw_free(in); fftw_free(out);
    return 0;
}

