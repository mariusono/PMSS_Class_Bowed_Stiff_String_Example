#include "MainComponent.h"

// OTHER FUNCTIONS

double linearMapping(float rangeIn_top, float rangeIn_bottom, float rangeOut_top, float rangeOut_bottom, float value) {
    double newValue = rangeOut_bottom + ((rangeOut_top - rangeOut_bottom) * (value - rangeIn_bottom) / (rangeIn_top - rangeIn_bottom));
    return newValue;
}


double exponentialMapping(float rangeIn_top, float rangeIn_bottom, float rangeOut_top, float rangeOut_bottom, float fac, float value)
{
    // make sure values passed to function are within the rangeIn_bottom rangeIn_top interval !!
    // maybe add an error exception here..
    // first map value to rangeIn to 0 - 1
    double valueMapped = 0.0 + ((1.0 - 0.0) * (value - rangeIn_bottom) / (rangeIn_top - rangeIn_bottom));
    
    // map to an exponential curve between 0 and 1 with a factor fac
    double mapToExp = (exp(valueMapped * fac) - 1) / (exp(fac)-1);

    // map back to desired output range
    double newValue = rangeOut_bottom + ((rangeOut_top - rangeOut_bottom) * (mapToExp - 0.0) / (1.0 - 0.0));

    return newValue;
}

double clamp(double in, double min, double max)
{
    if (in > max)
        return max;
    else if (in < min)
        return min;
    else
        return in;
}

std::vector<float> gOutputStrings;

//==============================================================================
MainComponent::MainComponent() : minOut(-1.0), maxOut(1.0)
{
    // Make sure you set the size of the component after
    // you add any child components.
    setSize((int) widthScreen, (int) heightScreen);

    // specify the number of input and output channels that we want to open
    setAudioChannels(0, 2);
    
    // remove mouse cursor from this
    setMouseCursor(MouseCursor::NoCursor);

    // SLIDER STUFF:

    addAndMakeVisible(dampingSlider);
    dampingSlider.setSliderStyle(juce::Slider::Rotary);
    dampingSlider.setTextBoxStyle(juce::Slider::NoTextBox, false, 0, 0); // Add a textbox..
    dampingSlider.setTextValueSuffix(" [-]");
    dampingSlider.setRange(0.1, 10.0);          // maybe scale these values afterwards
    dampingSlider.setValue(4.5); // with exp mapping this will result in 1.0ish
    dampingSlider.addListener(this);

    dampingLabel.setText("Freq Indep Damping", juce::dontSendNotification);
    dampingLabel.attachToComponent(&dampingSlider, true); // [4]
    addAndMakeVisible(dampingLabel);
    
    
    addAndMakeVisible(freqDampingSlider);
    freqDampingSlider.setSliderStyle(juce::Slider::Rotary);
    freqDampingSlider.setTextBoxStyle(juce::Slider::NoTextBox, false, 0, 0); // Add a textbox..
    freqDampingSlider.setTextValueSuffix(" [-]");
    freqDampingSlider.setRange(0.001, 0.15);          // maybe scale these values afterwards
    freqDampingSlider.setValue(0.085); // ADD EXP MAPPING ! with exp mapping this will result in 10.0ish
    freqDampingSlider.addListener(this);

    freqDampingLabel.setText("Freq Dep Damping", juce::dontSendNotification);
    freqDampingLabel.attachToComponent(&freqDampingSlider, true); // [4]
    addAndMakeVisible(freqDampingLabel);
    
    
    addAndMakeVisible(frParamSlider);
    frParamSlider.setSliderStyle(juce::Slider::Rotary);
    frParamSlider.setTextBoxStyle(juce::Slider::NoTextBox, false, 0, 0);
    frParamSlider.setTextValueSuffix(" [-]");
    frParamSlider.setRange(1.0, 15000.0);
    frParamSlider.setValue(7500.0); // with exp mapping this will result in 100.0ish
    frParamSlider.addListener(this);

    frParamLabel.setText("Friction Amount", juce::dontSendNotification);
    frParamLabel.attachToComponent(&frParamSlider, true); // [4]
    addAndMakeVisible(frParamLabel);
    
    
    addAndMakeVisible(volumeSlider);
    volumeSlider.setSliderStyle(juce::Slider::Rotary);
    volumeSlider.setTextBoxStyle(juce::Slider::NoTextBox, false, 0, 0);
    volumeSlider.setTextValueSuffix(" [-]");
    volumeSlider.setRange(-10.0, 10.0); // in dB. conversion done in sliderValueChanged routine
    volumeSlider.setValue(-1.0);
    volumeSlider.addListener(this);

    volumeLabel.setText("Volume", juce::dontSendNotification);
    volumeLabel.attachToComponent(&volumeSlider, true);
    addAndMakeVisible(volumeLabel);
}

MainComponent::~MainComponent()
{
    // This shuts down the audio device and clears the audio source.
    shutdownAudio();
}


void MainComponent::sliderValueChanged(Slider* slider)
{
    if (slider == &dampingSlider)
    {
        gSig0 = dampingSlider.getValue();
        gSig0 = exponentialMapping(10.0, 0.1, 10.0, 0.1, 4, gSig0);
        
        Logger::getCurrentLogger()->outputDebugString("gSig0: (" + String(gSig0) + ")");
    }
    else if (slider == &freqDampingSlider)
    {
        gSig1 = freqDampingSlider.getValue();
        gSig1 = exponentialMapping(0.15, 0.001, 0.15, 0.001, 2, gSig1);

////        gFrParam = exponentialMapping(1500.0, 1.0, 1500.0, 1.0, 4.0, gFrParam);
//        gFrParam = exponentialMapping(15000.0, 1.0, 15000.0, 1.0, 8.0, gFrParam);
//        gSig1 = linearMapping(15000.0, 1.0, 2.0, 0.0, gFrParam);

//        gSig1 = linearMapping(15000.0, 1.0, 1, 15000.0, gFrParam);
//        gSig1 = exponentialMapping(1.0, 15000.0, 1.0, 15000.0, -8.0, gFrParam);

        Logger::getCurrentLogger()->outputDebugString("gSig1: (" + String(gSig1) + ")");
    }
    else if (slider == &frParamSlider)
    {
        gFrParam = frParamSlider.getValue();
////        gFrParam = exponentialMapping(1500.0, 1.0, 1500.0, 1.0, 4.0, gFrParam);
//        gFrParam = exponentialMapping(15000.0, 1.0, 15000.0, 1.0, 8.0, gFrParam);
        gStickFact = linearMapping(15000.0, 1.0, 2.0, 0.0, gFrParam);

        gFrParam = linearMapping(15000.0, 1.0, 1, 15000.0, gFrParam);
        gFrParam = exponentialMapping(1.0, 15000.0, 1.0, 15000.0, -8.0, gFrParam);

        Logger::getCurrentLogger()->outputDebugString("gFrParam: (" + String(gFrParam) + ")");
    }
    else if (slider == &volumeSlider)
    {
        gVolume = volumeSlider.getValue();
        gVolume = powf(10.0, gVolume / 20);
        if (gVolume < 0.317)
        {
            gVolume = 0;
        }
        
        Logger::getCurrentLogger()->outputDebugString("gVolume: (" + String(gVolume) + ")");

    }
}


//==============================================================================
void MainComponent::prepareToPlay (int samplesPerBlockExpected, double sampleRate)
{
    
    // ADD SENSEL
    for (int i = 0; i < amountOfSensels; ++i)
    {
        sensels.add(new Sensel(i));
        Logger::getCurrentLogger()->outputDebugString("Sensel added");
    }

    hiResCallbackFreq = 150.0;
    // start the hi-res timer
    if (sensels.size() != 0)
    {
        if (sensels[0]->senselDetected)
        {
            HighResolutionTimer::startTimer(1000.0 / hiResCallbackFreq); // 150 Hz
        }
    }
    
    globalCurrentSample = 0;

    fs = sampleRate;
    bufferSize = samplesPerBlockExpected;

    k = 1/fs; // time step
    
    L = 0.5;
    rho = 7850;
    r = 5e-4;
    T = 1000;
    E = 2e11;
    sig0 = 1;
//    sig1 = 0.01;
    sig1 = 0.151; // max possible val

    A = double_Pi*r*r; // [m^2]
    I = double_Pi*r*r*r*r/4;
    K = sqrt(E*I/(rho*A));
    

    // Bowing params
    
    FB = 0.2; // [N]
    vB = 0.2; // bow velocity (m/s)
    a = 80.0; // friction law free parameter (1/m^2) (a) % decreasing sig increases the stick time
    tol = 1e-7; // tolerance for Newton-Raphson method
    A_NR = sqrt(2*a)*exp(0.5);
    
    inp_bow_pos_x = 0.5;
    
    stickFact = 1;
    
    fN_var = 0.0;
    vB_var = 0.0;
    
//    fN_var = 0.2;
//    vB_var = 0.1;
    
    freq = 110.0;
//    freq = 440.0;

    // Derived params for num sim
    c = freq*2*L; // wave speed term
    h = sqrt((c * c * k * k + 4 * sig1 * k + sqrt(pow((c * c * k * k + 4 * sig1 * k),2) + 16 * K * K * k * k)) / 2);
    N = floor(L / h);
    h = L / N;
    
    // Bow interpolation and spreading function init:
    I_B.resize(N, 0);
    J_B.resize(N, 0);

    // INITIALIZE STRING STATE VECTORS
    uVecs.reserve(3);

    for (int i = 0; i < 3; ++i)
        uVecs.push_back(std::vector<double>(N, 0));

    u.resize(3);

    for (int i = 0; i < u.size(); ++i)
        u[i] = &uVecs[i][0];
    
    
    
    preOutputVec.resize(1, 0);
    gOutputStrings.resize(1, 0);
    keyDownVec.resize(1, false);

//    activeStrings.resize(polyphony, nullptr);

    // Timer for graphics
    startTimerHz(72);

    gPosition.resize(2, 0);
    gPositionRaw.resize(2, 0);
    gPositionPrev.resize(2, 0);
    
    
    // position of bow ?
    x_inp_var = 0.3; // bowing pos in percentage
    y_inp_var = 0.5;
    

}

void MainComponent::getNextAudioBlock (const juce::AudioSourceChannelInfo& bufferToFill)
{
    sig0 = gSig0;
    sig1 = gSig1;
    a = gFrParam;
    stickFact = gStickFact;
//    fN_var = 0.1;
//    vB_var = 0.2;

    
    // update volume at every buffersize to avoid clicks.. (maybe)
    double volPerBuffer = gVolume;
    //    Logger::getCurrentLogger()->outputDebugString("output: (" + String(gOutputStrings[0]) + ")");
    
//    Logger::getCurrentLogger()->outputDebugString("output: (" + String(gOutputStrings[0]) + ")");
    
    
    // Your audio-processing code goes here!
    for (int channel = 0; channel < bufferToFill.buffer->getNumChannels(); ++channel)
    {
        float *const channelData = bufferToFill.buffer->getWritePointer(channel, bufferToFill.startSample);

        if (channel == 0)
        {


            for (int i = 0; i < bufferToFill.buffer->getNumSamples(); i++)
            {
                bp = x_inp_var;
                bP = floor(bp * N - 1); // -1 for alignment with Matlab..
                alpha_bow = bp * N - 1 - bP;

                int bP_m1 = bP - 1;
                int bP_p1 = bP + 1;
                int bP_p2 = bP + 2;

                std::fill(I_B.begin(), I_B.end(), 0); // fill with zeros

                I_B[bP_m1] = (alpha_bow * (alpha_bow - 1) * (alpha_bow - 2)) / -6.0;
                I_B[bP] = ((alpha_bow - 1) * (alpha_bow + 1) * (alpha_bow - 2)) / 2.0;
                I_B[bP_p1] = (alpha_bow * (alpha_bow + 1) * (alpha_bow - 2)) / -2.0;
                I_B[bP_p2] = (alpha_bow * (alpha_bow + 1) * (alpha_bow - 1)) / 6.0;

                for (int i = 0; i < I_B.size(); ++i)
                {
                    J_B[i] = I_B[i] * (1 / h); // speed up: keep divisions out of loop !
                }
                
                double I_B_J_B = 0;
                double I_B_u = 0;
                double I_B_uPrev = 0;
                double I_B_dxx_u = 0;
                double I_B_dxx_uPrev = 0;
                double I_B_dxxxx_u = 0;
                int idx_p1;
                int idx_p2;
                int idx_m1;
                int idx_m2;

                for (int idx = 2; idx < N-2; ++idx)
                {
                    idx_p1 = idx + 1;
                    idx_p2 = idx + 2;
                    idx_m1 = idx - 1;
                    idx_m2 = idx - 2;

                    I_B_J_B = I_B_J_B + I_B[idx] * J_B[idx];
                    I_B_u = I_B_u + I_B[idx] * u[1][idx];
                    I_B_uPrev = I_B_uPrev + I_B[idx] * u[2][idx];

                    I_B_dxx_u = I_B_dxx_u + I_B[idx] * (u[1][idx_p1] - 2. * u[1][idx] + u[1][idx_m1]) * (1 / (h * h));
                    I_B_dxx_uPrev = I_B_dxx_uPrev + I_B[idx] * (u[2][idx_p1] - 2. * u[2][idx] + u[2][idx_m1]) * (1 / (h * h));

                    I_B_dxxxx_u = I_B_dxxxx_u
                        + I_B[idx] * (u[1][idx_p2] - 4. * u[1][idx_p1] + 6. * u[1][idx] - 4. * u[1][idx_m1] + u[1][idx_m2]) * (1 / (h * h * h * h));
                }
                
                double q = (-2 / k) * (1 / k) * (I_B_u - I_B_uPrev)
                            + (2 / k) * vB_var
                            + 2 * sig0 * vB_var
                            - (c * c) * I_B_dxx_u
                            + K * K * I_B_dxxxx_u
                            - (2 * sig1) * (1 / k) * I_B_dxx_u
                            + (2 * sig1) * (1 / k) * I_B_dxx_uPrev;
                
                if (fN_var > 0)
                {
                    double eps = 1;
                    //w_rnd_last = -1 + (1 - (-1)).*rand(1);

                    int iter_check = 0;
                    vrel_last = -vB; // should this be vB_var ?
                    // Newton-Raphson iterative scheme
                    while ((eps > tol))
                    {
                        ++iter_check;

                        g1 = (2/k + 2*sig0) * vrel_last + stickFact*I_B_J_B * fN_var * sqrt(2 * a) * exp(0.5) * vrel_last * exp(- a *vrel_last * vrel_last) + q;

                        g1_d_vrel = stickFact*I_B_J_B * (fN_var * sqrt(2 * a) * exp(0.5) * exp(-a * vrel_last * vrel_last) * (1 + 2 * a * vrel_last * vrel_last)) / (rho * A) + 2 * sig0 + 2 / k;
                        
                        vrel = vrel_last - g1 / g1_d_vrel;
                        eps = abs(vrel - vrel_last);

                        vrel_last = vrel;
                        if (iter_check == 99)
                            //Logger::getCurrentLogger()->outputDebugString("iter_check: (" + String(iter_check) + ")");
                            break;

                    }
                    vrel = vrel_last;

                    Fbow = stickFact*fN_var * exp(0.5) * sqrt(2 * a) * vrel *exp(-a * vrel * vrel);

                    //Fbow = fN_var * A_fr_model * vrel * exp(-sig * vrel * vrel);
                    phi_vrel = exp(0.5) * sqrt(2 * a) * vrel * exp(-a * vrel * vrel);
                }
                else
                {
                    Fbow = 0;
                    phi_vrel = exp(0.5) * sqrt(2 * a) * vrel * exp(-a * vrel * vrel);
                    vrel_last = -vB;
                }
                
                // Update equations
                for (int idx = 2; idx < N - 2; ++idx)
                {
                    idx_p1 = idx + 1;
                    idx_p2 = idx + 2;
                    idx_m1 = idx - 1;
                    idx_m2 = idx - 2;

                     // maybe you can reused the definitions from above I_S_etc

                    u[0][idx] = (k * k / (1 + sig0 * k)) * ((c * c / (h * h)) * (u[1][idx_p1] - 2 * u[1][idx] + u[1][idx_m1])
                            + (2 * sig1 / (k * h * h)) * (u[1][idx_p1] - 2 * u[1][idx] + u[1][idx_m1] - u[2][idx_p1] + 2 * u[2][idx] - u[2][idx_m1])
                            - (K * K) * (1 / (h * h * h * h)) * (u[1][idx_p2] - 4 * u[1][idx_p1] + 6 * u[1][idx] - 4 * u[1][idx_m1] + u[1][idx_m2])
                            - J_B[idx] * Fbow / (rho * A)
                            + (2 / (k * k)) * u[1][idx] - (1 - sig0 * k) * u[2][idx] / (k * k));
                }
                
                u[0][1] = (k * k / (1 + sig0 * k)) * ((c * c / (h * h)) * (u[1][2] - 2 * u[1][1] + u[1][0])
                    + (2 * sig1 / (k * h * h)) * (u[1][2] - 2 * u[1][1] + u[1][0] - u[2][2] + 2 * u[2][1] - u[2][0])
                    - (K * K) * (1 / (h * h * h * h)) * (u[1][3] - 4 * u[1][2] + 6 * u[1][1] - 4 * u[1][0] - u[1][1])
                    - J_B[1] * Fbow / (rho * A)
                    + (2 / (k * k)) * u[1][1] - (1 - sig0 * k) * u[2][1] / (k * k));

                int idx = N - 2;
                idx_p1 = idx + 1;
                idx_p2 = idx + 2;
                idx_m1 = idx - 1;
                idx_m2 = idx - 2;

                u[0][idx] = (k * k / (1 + sig0 * k)) * ((c * c / (h * h)) * (u[1][idx_p1] - 2 * u[1][idx] + u[1][idx_m1])
                    + (2 * sig1 / (k * h * h)) * (u[1][idx_p1] - 2 * u[1][idx] + u[1][idx_m1] - u[2][idx_p1] + 2 * u[2][idx] - u[2][idx_m1])
                    - (K * K) * (1 / (h * h * h * h)) * (-u[1][idx] - 4 * u[1][idx_p1] + 6 * u[1][idx] - 4 * u[1][idx_m1] + u[1][idx_m2])
                    - J_B[idx] * Fbow / (rho * A)
                    + (2 / (k * k)) * u[1][idx] - (1 - sig0 * k) * u[2][idx] / (k * k));

                int idx_out = round(N * 0.7);

                double vel_out = (1/(2*k)) * (u[2][idx_out] - u[0][idx_out]);
                
                uTmp = u[2];
                u[2] = u[1];
                u[1] = u[0];
                u[0] = uTmp;
                
//                                  double output = u[0][idx_out];
                double output = vel_out;



                
                float gainVal = 0.1;
                output = output * gainVal;
                output = output * volPerBuffer; // add global volume contribution
            
                float stringSound = output;
        
                gOutputStrings[0] = stringSound;
                
                if (output > maxOut) // output is limited..
                {
                    Logger::getCurrentLogger()->outputDebugString("Output is too large!");
                    output = maxOut;
                }
                else if (output < minOut) {
                    Logger::getCurrentLogger()->outputDebugString("Output is too small!");
                    output = minOut;
                }
                channelData[i] = output;
            }
        }
        else
        {
            memcpy(channelData,
                bufferToFill.buffer->getReadPointer(0),
                bufferToFill.numSamples * sizeof(float));
        }
    }

//    bufferToFill.clearActiveBufferRegion();
}

//=============================================================================
void MainComponent::releaseResources()
{
    // This will be called when the audio device stops, or when it is being
    // restarted due to a setting change.

    // For more details, see the help for AudioProcessor::releaseResources()
}

//==============================================================================
void MainComponent::paint(juce::Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    // You can add your drawing code here!

    Rectangle<int> area = getLocalBounds();

//    g.fillAll(getLookAndFeel().findColour(ResizableWindow::backgroundColourId));

    if (flagMouseUp)
    {
        g.setColour(Colours::grey);
        g.setOpacity(opa_level);
    }
    else
    {
        g.setColour(Colours::orange);
        g.setOpacity(opa_level);
    }

    int heightRect = 15;
    int xRect = gPositionRaw[0] - (area.getWidth()/2)/2;
    int yRect = gPositionRaw[1] - heightRect/2;
    int widthRect = area.getWidth()/2;
//    int heightRect = 30;

    g.fillRect(xRect,yRect,widthRect,heightRect);

    g.setOpacity(1.0);
    g.setColour(Colours::cyan);

    // draw the state
    g.strokePath(visualiseState (g, 400000), PathStrokeType(2.0f));

}


Path MainComponent::visualiseState (Graphics& g, double visualScaling)
{
    // String-boundaries are in the horizontal middle of the component
//    double stringBoundaries= getHeight() / 2.0;
    double stringBoundaries_y  = 0;

    // initialise path
    Path stringPath;

    // start path
    stringPath.startNewSubPath (stringBoundaries_y,u[1][0] * visualScaling);

    double spacing = getHeight() / static_cast<double>(N);
//    double y = spacing;
    double y = 0;

    for (int l = 0; l < N; l++) // if you don't save the boundaries use l < N
    {
        // Needs to be -u, because a positive u would visually go down
        float newX = u[1][l] * visualScaling + getWidth()/2;

        // if we get NAN values, make sure that we don't get an exception
        if (isnan(newX))
            newX = 0;

        stringPath.lineTo (newX,y);
        y += spacing;
    }
    // if you don't save the boundaries, and add a stringPath.lineTo (x, getWidth()) here to end the statedrawing

    return stringPath;
}


void MainComponent::resized()
{
    
    
    // This is called when the MainContentComponent is resized.
    // If you add any child components, this is where you should
    // update their positions.

    Rectangle<int> area = getLocalBounds();

    int x_dampSlider = area.getX() + area.getWidth() / 8;
    int y_dampSlider = area.getY() + area.getHeight() / 8;
    int width_dampSlider = 100;
    int height_dampSlider = 100;

    dampingSlider.setBounds(x_dampSlider,y_dampSlider,width_dampSlider,height_dampSlider);

    int x_frParamSlider = area.getX() + area.getWidth() / 8;
    int y_frParamSlider = 240 + area.getY() + area.getHeight() / 8;
    int width_frParamSlider = 100;
    int height_frParamSlider = 100;

    frParamSlider.setBounds(x_frParamSlider, y_frParamSlider, width_frParamSlider, height_frParamSlider);


    int x_volumeSlider = area.getX() + area.getWidth() / 8;
    int y_volumeSlider = 360 + area.getY() + area.getHeight() / 8;
    int width_volumeSlider = 100;
    int height_volumeSlider = 100;

    volumeSlider.setBounds(x_volumeSlider, y_volumeSlider, width_volumeSlider, height_volumeSlider);
    

    int x_freqDampSlider = area.getX() + area.getWidth() / 8;
    int y_freqDampSlider  = 120 + area.getY() + area.getHeight() / 8;
    int width_freqDampSlider = 100;
    int height_freqDampSlider = 100;

    freqDampingSlider.setBounds(x_freqDampSlider, y_freqDampSlider, width_freqDampSlider, height_freqDampSlider);
}



void MainComponent::timerCallback()
{
    repaint(); // update the graphics X times a second
}



void MainComponent::hiResTimerCallback()
{
//    repaint(); // update the graphics X times a second
    
    
    for (auto sensel : sensels)
    {
        if (sensel->senselDetected)
        {
            sensel->check();
            unsigned int fingerCount = sensel->contactAmount;
//            int index = sensel->senselIndex;
//            fingerCount = 1; // force to always be 0..
                        
            if (fingerCount > 1) // limit to only one possible touch..
            {
//                fingerCount = 1;
                flagBow_Sensel = true;
                flagMouseUp = false;
            }
            else if (fingerCount <= 1)
            {
                flagBow_Sensel = false;
                flagMouseUp = true;
            }
//
//            if (fingerCount == 1)
//            {
                float x = sensel->fingers[0].x;
                float y = sensel->fingers[0].y;
                
                x = linearMapping(0.0, 0.95, 0.0, 1.0, x);
                y = linearMapping(0.0, 0.93, 0.0, 1.0, y);
                x = clamp(x, 0.0, 1.0);
                y = clamp(y, 0.0, 1.0);

                gPosition[0] = x;
                gPosition[1] = y;
            
                x_inp_var = y;
                
                x_inp_var = clamp(x_inp_var, 0.1, 0.9);
            
                gPositionRaw[0] = linearMapping(0.0, 1.0, 0.0, widthScreen, x);
                gPositionRaw[1] = linearMapping(0.0, 1.0, 0.0, heightScreen, y);

                if (flagBow_Sensel)
                {
//                    ViolinStrings[j]->setBow(true);
//                    ViolinStrings[j]->activate();

                    isBowing = true;
                    active = true;
                }
                double VbRaw = (gPosition[0] - gPositionPrev[0]) * hiResCallbackFreq;
                
//                Logger::getCurrentLogger()->outputDebugString("Vb Raw: (" + String(VbRaw) + ")");

                
                double maxVb = 0.3;
                double Vb = linearMapping(-12.0, 12.0, -maxVb, maxVb, VbRaw);
                if (Vb > -0.03 && Vb < 0.03) // limit Vb to avoid small noise during small displacements
                {
                    Vb = 0.0;
                }
                if (Vb > 0.4) // limit Vb to avoid small noise during small displacements
                {
                    Vb = 0.4;
                }

//                Logger::getCurrentLogger()->outputDebugString("Vb: (" + String(Vb) + ")");
                //Logger::getCurrentLogger()->outputDebugString("FB: (" + String(mass_spring_dampers[idx]->getFb()) + ")");

                vB_var = Vb;
                Logger::getCurrentLogger()->outputDebugString("vB_var: (" + String(vB_var) + ")");
//                ViolinStrings[0]->setVb(Vb);
//                ViolinStrings[0]->setVb(0.2);

                double forceRead = sensel->fingers[0].force;
//                    Logger::getCurrentLogger()->outputDebugString("forceRead: (" + String(forceRead) + ")");

                forceRead = clamp(forceRead, 0.0, 0.2);
                
                double maxFB = 1.0;
//                double FB = 800;
                double FB = linearMapping(0.0, 0.15, 0.0, maxFB, forceRead);
                FB = clamp(FB, 0.0, maxFB);
                
                if (vB_var == 0)
                {
                    fN_var = 0;
                }
                else
                {
                    fN_var = FB;
                }

//                Logger::getCurrentLogger()->outputDebugString("FB: (" + String(FB) + ")");

            //    mass_spring_dampers[idx_mass]->setFb(FB); // update force on each mass..
                
//                for (int j = 0; j < numStrings; ++j) // update force on all masses at once (better!)
//                {
//                    ViolinStrings[j]->setFb(FB);
////                    ViolinStrings[j]->setFb(0.1);
//                }
            
//                Logger::getCurrentLogger()->outputDebugString("FB: (" + String(ViolinStrings[0]->getFb()) + ")");
//                Logger::getCurrentLogger()->outputDebugString("vB: (" + String(ViolinStrings[0]->getVb()) + ")");
                

            
                opa_level = FB/maxFB;
                if (opa_level < 0.1)
                {
                    opa_level = 0.1;
                }
                else if (opa_level > 0.9)
                {
                    opa_level = 0.9;
                }
                
                gPositionPrev[0] = gPosition[0];
                gPositionPrev[1] = gPosition[1];
            
            if (!flagBow_Sensel)
            {
                gPositionPrev[0] = gPosition[0];
                gPositionPrev[1] = gPosition[1];
                
//                    ViolinStrings[j]->setBow (false);
                isBowing = false;
                fN_var = 0;
                vB_var = 0;
            }
            if (fingerCount == 0)
            {
                flagMouseUp = true;
            }
        }
    }

}
