#pragma once

#include <JuceHeader.h>
#include "SenselWrapper.h"

//==============================================================================
/*
    This component lives inside our window, and this is where you should put all
    your controls and content.
*/
class MainComponent  : public juce::AudioAppComponent,
                        public juce::HighResolutionTimer,
                        public juce::Timer,
                        public juce::Slider::Listener,
                        public juce::MouseListener
{
public:
    //==============================================================================
    MainComponent();
    ~MainComponent() override;

    //==============================================================================
    void prepareToPlay (int samplesPerBlockExpected, double sampleRate) override;
    void getNextAudioBlock (const juce::AudioSourceChannelInfo& bufferToFill) override;
    void releaseResources() override;

    //==============================================================================
    void paint (juce::Graphics& g) override;
    void resized() override;
    
    // function to draw the state of the string
    Path visualiseState (Graphics& g, double visualScaling);

    // Hi resolution timer callback
    void hiResTimerCallback() override;

    // Graphics timer callback
    void timerCallback() override;

    int idx_mass = 0;
    double wheelDeltaY = 0.0;

    // Slider stuff
    void sliderValueChanged(Slider* slider) override;

private:
    //==============================================================================
    // Your private member variables go here...

    
    
    // init params
    double fs, freq, mass, sig0, sig1, vB, A_NR, tol, a, k, h, c;
    double FB;
    double L,rho,r,T,E;
    double A,I,K;
    double inp_bow_pos_x; // in perc
    double N;
    
    double fN_var, vB_var;
    double Fbow;
    double phi_vrel;
    
    double stickFact;
    
    // BOWING POSITION INIT // this might go in MainComponent..
    /// Parameters you may want to modulate:
    double bp = 0.25; // in percentage
    int bP;
    double alpha_bow;
    
    // Bowing interpolation fct and spreading fct
    std::vector<double> I_B; // interpolant grid for bowing pos
    std::vector<double> J_B; // spreading function grid for bowing pos
    
    // pointers to STRING states
    std::vector<double*> u;
    
    // Vector states
    std::vector<std::vector<double>> uVecs;
    double* uTmp = nullptr;
    
    
    // INITIALIZE PARAMS FOR NEWTON-RAPHSON
    double vrel;
    double vrel_last;
    double g1;
    double g1_d_vrel;
    
    
    // process params
    double F_fr;
    
    // NR params
    double eps, numerator, denominator;
    
    // State params
    double u_n, u_nm1, u_np1; // velocity and displacement of mass

    
//    gamma, k, s0, s1, B, kappa, h, N, lambdaSq, muSq, kOh, gOh, a, BM, Vb, Fb, pickup, tol, q, qPrev, bp, b, eps;
    bool isBowing = false;
//    std::vector<double> u, uPrev, uNext;
    bool active = false;
    
    unsigned long countGlobal;
//    unsigned long t = 0;
    
    
    double bufferSize;
    float minOut;
    float maxOut;
    
    int testCommit;
    double globalCurrentSample;

    //std::vector<std::unique_ptr<mass_spring_damper>> mass_spring_dampers;
    //std::vector<std::unique_ptr<mass_spring_damper>> activeMasses;

    int numStrings;
    int octave;
    double FB_Max;
    

    double x_inp_var;
    double y_inp_var;
    
    
    //std::vector<const char> keys = {'A', 'W', 'S', 'E', 'D', 'F', 'T', 'G', 'Y', 'H', 'U', 'J', 'K'}; // const char does not work in this version of VS ... !!!
//    std::vector<char> keys = { 'A', 'W', 'S', 'E', 'D', 'F', 'T', 'G', 'Y', 'H', 'U', 'J', 'K' };

    int polyphony;
    int currentPoly = 0;

    std::vector<double> preOutputVec;
    std::vector<bool> keyDownVec;

    double eYPrev = 0;
    long int timePrev = 0;;
    
    std::vector<double> gPosition;
    std::vector<double> gPositionRaw;
    std::vector<double> gPositionPrev;

    // SLIDER STUFF
    Slider dampingSlider;
    Label dampingLabel;

    Slider freqDampingSlider;
    Label freqDampingLabel;

    
    Slider frParamSlider;
    Label frParamLabel;
    
    Slider volumeSlider;
    Label volumeLabel;
    
    double gDampingVal;
    double gFrParam;
    double gSig0;
    double gSig1;
    double gVolume;
    double gStickFact;

    // For graphics
    // Opacity
    float opa_level = 1.0;
    float widthScreen = 800;
    float heightScreen = 600;
    bool flagMouseUp = true;

    // Sensel Stuff
    OwnedArray<Sensel> sensels;
    int amountOfSensels = 1;
    
    double hiResCallbackFreq;

    double gIdx_map;
    bool flagBow_Sensel;
    
    
    // pointers to STRING states
    std::vector<double*> u_sel;
    
    
    
    // FOR DC BLOCKER
    double y_nm1;
    double y_n;
    double x_nm1;
    double x_n;
    double R_fac;
    
    // FOR LIMITER
    double at;
    double rt;
    int delaySamples;
    double outPeak;
    double g_lim;
    double coeff;
    double lt;
    double outputLim;
    int BUFFER_SIZE = 1024;
    
    int gReadPointer_limiter;
    int gWritePointer_limiter;
    std::vector<double> delayLine_limiter;


    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainComponent)
};
