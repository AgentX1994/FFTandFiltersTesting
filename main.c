#include <complex.h>
#include <fftw3.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <sndfile.h>
#include <portaudio.h>

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

// Struct to hold return values from chebyshev_poles function
typedef struct {
    double a0;
    double a1;
    double a2;
    double b1;
    double b2;
} POLES;

// Calculates the poles struct for chebyshev function
// double fc:   filter cutoff (0 to 0.5 percent of sampling frequency)
// bool lp:     false for low pass, true for high pass
// double pr:   Percent Ripple (0 to 29)
// int np:      Number of poles (2 to 20, must be even)
// int p:       Current pole being calculated
POLES chebyshev_poles(double fc, bool lp, double pr, int np, int p)
{
    // Calculate the location of the pole on the unit circle
    double angle = M_PI/(np*2) + (p-1)*M_PI/np;
    double rp = - cos(angle);
    double ip =   sin(angle);

    double es = 0;
    double vx = 0;
    double kx = 0;
    if (pr != 0.0) // Warp from circle to eclipse
    {
        double temp = 100/(100 - pr);
        es = sqrt(temp*temp - 1);
        vx = log((1/es) + sqrt(1/(es*es) + 1))/np;

        kx = log((1/es) + sqrt(1/(es*es) - 1))/np;
        kx = (exp(kx) + exp(-kx)) / 2;
        
        rp *= (exp(vx) - exp(-vx))/(2*kx);
        ip *= (exp(vx) + exp(-vx))/(2*kx);
    } 

    // s domain to z domain transform
    double t = 2*tan(1.0/2.0);
    double w = 2*M_PI*fc;
    double m = rp*rp + ip*ip;
    double d = 4 - 4*rp*t + m*t*t;

    double x0 = t*t/d;
    double x1 = 2*x0;
    double x2 = x0;
    double y1 = (8 - 2*m*t*t)/d;
    double y2 = (-4 - 4*rp*t - m*t*t)/d;

    // LP to LP or LP to HP transform
    double k;
    if (lp)
    {
        k = - cos(w/2 + 1.0/2.0) / cos(w/2 - 1.0/2.0);
    }
    else
    {
        k = sin(1.0/2.0 - w/2) / sin(1.0/2.0 + w/2);
    }

    d = 1 + y1*k - y2*k*k;
    
    POLES ret;
    ret.a0 = (x0 - x1*k + x2*k*k)/d;
    ret.a1 = (-2*x0*k + x1 + x1*k*k - 2*x2*k)/d;
    ret.a2 = (x0*k*k - x1*k + x2)/d;
    ret.b1 = (2*k + y1 + y1*k*k - 2*y2*k)/d;
    ret.b2 = (-k*k - y1*k + y2)/d;

    if (lp)
    {
        ret.a1 = -ret.a1;
        ret.b1 = -ret.b1;
    }

    return ret;
}

// Struct to hold return values from chebyshev (filter coefficients)
typedef struct {
    double *a;
    double *b;
} FLT_COEFF;

// Chebyshev filter coefficients
// parameters:
// double fc:   filter cutoff (0 to 0.5 percent of sampling frequency)
// bool lp:     false for low pass, true for high pass
// double pr:   Percent Ripple (0 to 29)
// int np:      Number of poles (2 to 20, must be even)
// double *a:   output buffer to hold a coefficients
// double *b:   output buffer to hold b coefficients
FLT_COEFF chebyshev(float fc, bool lp, double pr, int np)
{
    int n = 23; // size of all arrays
    FLT_COEFF c;
    c.a = malloc(sizeof(double)*n);
    c.b = malloc(sizeof(double)*n);
    double ta[n];
    double tb[n];
    for (int i = 0; i < n; ++i)
    {
        c.a[i] = 0;
        c.b[i] = 0;
    }
    c.a[2] = 1;
    c.b[2] = 1;

    // Loop for each pole pair
    for (int p = 1; p < np/2+1; ++p)
    {
        POLES ret = chebyshev_poles(fc, lp, pr, np, p);

        // Add coefficients to the cascade
        for (int i = 0; i < n; ++i)
        {
            ta[i] = c.a[i];
            tb[i] = c.b[i];
        }

        for (int i = 2; i < n; ++i)
        {
            c.a[i] = ret.a0*ta[i] + ret.a1*ta[i-1] + ret.a2*ta[i-2];
            c.b[i] =        tb[i] - ret.b1*tb[i-1] - ret.b2*tb[i-2];
        }

    }

    // Finish combining coefficients
    c.b[2] = 0;
    for (int i = 0; i < n-2; ++i)
    {
        c.a[i] =  c.a[i+2];
        c.b[i] = -c.b[i+2];
    }

    // Normalize the gain
    double sa = 0;
    double sb = 0;
    for (int i = 0; i < n; ++i)
    {
        if (!lp)
        {
            sa += c.a[i];
            sb += c.b[i];
        }
        else
        {
            sa += c.a[i]*pow(-1,i); 
            sb += c.b[i]*pow(-1,i);
        }
    }
    double gain = sa/(1-sb);
    for (int i = 0; i < n; ++i)
        c.a[i] /= gain;
    
    return c;
}

// Applies filter to input buffer and stores in output buffer
// Calculates the convolution of in by a, along with rescursive feedback from b
// params:
// double *in               input buffer
// int in_n                 input buffer length
// double *a                filter kernel
// double *b                filter kernel (recursive part)
// int flt_n                filter kernel length (assumes equal length for a and b)
// double *out              output buffer
// int out_n                output buffer length
void apply_filter(double *in, int in_n, double *a, double *b, int flt_n, double *out, int out_n)
{
    for (int i = 0; i < out_n; ++i)
    {
        out[i] = 0;
        for (int j = 0; j < flt_n; ++j)
        {
            double atemp = 0;
            double btemp = 0;
            if (i - j >= 0)
            {
                btemp = b[j]*out[i-j];
                if (i - j < in_n)
                {
                    atemp = a[j]*in[i-j];
                }
            }
            out[i] += atemp + btemp;
        }
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

