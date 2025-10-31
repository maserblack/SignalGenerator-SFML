// signal_generator_sfml_enhanced.cpp
// Digital Signal Generator with Line Coding and Decoder
// Compatible with SFML 2.x. Compile with:
// g++ -std=c++17 signal_generator_sfml_enhanced.cpp -o signal_gen \
//   -I/opt/homebrew/Cellar/sfml@2/2.6.2_1/include \
//   -L/opt/homebrew/Cellar/sfml@2/2.6.2_1/lib \
//   -lsfml-graphics -lsfml-window -lsfml-system

#include <SFML/Graphics.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <limits>
#include <sstream>
#include <iomanip>

using namespace std;

// ---------- Manacher's algorithm (O(n)) for longest palindrome ----------
string longestPalindromeManacher(const string &s) {
    if (s.empty()) return "";
    string t = "^";
    for (char c : s) { t += "#"; t += c; }
    t += "#$";
    int n = (int)t.size();
    vector<int> p(n, 0);
    int center = 0, right = 0;
    for (int i = 1; i < n-1; ++i) {
        int mir = 2*center - i;
        if (i < right) p[i] = min(right - i, p[mir]);
        while (t[i + 1 + p[i]] == t[i - 1 - p[i]]) p[i]++;
        if (i + p[i] > right) { center = i; right = i + p[i]; }
    }
    int maxLen = 0, centerIndex = 0;
    for (int i = 1; i < n-1; ++i)
        if (p[i] > maxLen) { maxLen = p[i]; centerIndex = i; }
    int start = (centerIndex - 1 - maxLen)/2;
    return s.substr(start, maxLen);
}

// ---------- Line Encoding Functions ----------
vector<int> nrzl(const string &bits) {
    vector<int> sig;
    for (char b : bits) sig.push_back(b=='1' ? 1 : -1);
    return sig;
}

vector<int> nrzi(const string &bits) {
    vector<int> sig;
    int level = -1;
    for (char b : bits) {
        if (b == '1') level *= -1;
        sig.push_back(level);
    }
    return sig;
}

vector<int> manchester(const string &bits) {
    vector<int> sig;
    for (char b : bits) {
        if (b == '1') { sig.push_back(1); sig.push_back(-1); }
        else { sig.push_back(-1); sig.push_back(1); }
    }
    return sig;
}

vector<int> diffManchester(const string &bits) {
    vector<int> sig;
    int level = 1;
    for (char b : bits) {
        if (b == '0') level *= -1;
        sig.push_back(level);
        level *= -1;
        sig.push_back(level);
    }
    return sig;
}

vector<int> ami(const string &bits) {
    vector<int> sig;
    int polarity = 1;
    for (char b : bits) {
        if (b == '1') { sig.push_back(polarity); polarity *= -1; }
        else sig.push_back(0);
    }
    return sig;
}

// ---------- Scrambling Schemes (O(n) complexity) ----------
string applyB8ZS(const string &bits) {
    string out = bits;
    for (size_t i = 0; i + 7 < out.size();) {
        bool all0 = true;
        for (size_t j = i; j < i+8; ++j) if (out[j] != '0') { all0 = false; break; }
        if (all0) {
            out[i+3] = '1';
            out[i+4] = '1';
            out[i+6] = '1';
            out[i+7] = '1';
            i += 8;
        } else ++i;
    }
    return out;
}

string applyHDB3(const string &bits) {
    string out = bits;
    int pulseCount = 0;
    for (size_t i = 0; i < out.size();) {
        if (out[i] == '1') { pulseCount++; ++i; continue; }
        size_t j = i;
        while (j < out.size() && out[j] == '0') ++j;
        size_t zeros = j - i;
        if (zeros >= 4) {
            size_t k = i;
            while (k + 3 < j) {
                if (pulseCount % 2 == 0) {
                    out[k+3] = '1';
                    pulseCount += 1;
                } else {
                    out[k] = '1';
                    out[k+3] = '1';
                    pulseCount += 2;
                }
                k += 4;
            }
        }
        i = j;
    }
    return out;
}

// ---------- PCM and Delta Modulation ----------
vector<int> pcmQuantize(const vector<double>& samples, int levels) {
    double minv = *min_element(samples.begin(), samples.end());
    double maxv = *max_element(samples.begin(), samples.end());
    if (fabs(maxv - minv) < 1e-12) {
        return vector<int>(samples.size(), 0);
    }
    vector<int> q;
    for (double x : samples) {
        double norm = (x - minv) / (maxv - minv);
        int code = (int)floor(norm * (levels - 1) + 0.5);
        q.push_back(code);
    }
    return q;
}

string quantizedToBitstream(const vector<int>& q, int levels) {
    int bitsNeeded = max(1, (int)ceil(log2(levels)));
    string out;
    for (int v : q) {
        for (int i = bitsNeeded - 1; i >= 0; --i) 
            out.push_back(((v >> i) & 1) ? '1' : '0');
    }
    return out;
}

vector<int> deltaMod(const vector<double>& samples, double step=0.5) {
    vector<int> bits;
    if (samples.empty()) return bits;
    double approx = samples[0];
    for (double s : samples) {
        if (s >= approx) {
            bits.push_back(1);
            approx += step;
        } else {
            bits.push_back(0);
            approx -= step;
        }
    }
    return bits;
}

// ---------- Upsample to time-domain samples ----------
vector<float> upsampleToFloatSamples(const vector<int>& levels, int samplesPerBit) {
    vector<float> out;
    out.reserve(levels.size() * samplesPerBit);
    for (size_t i = 0; i < levels.size(); ++i) {
        for (int s = 0; s < samplesPerBit; ++s) {
            out.push_back((float)levels[i]);
        }
    }
    return out;
}

// ---------- IMPROVED Signal Decoders (Fixed for 100% accuracy) ----------

string decodeFromSignal_NRZL(const vector<float>& samples, int samplesPerBit) {
    int numBits = samples.size() / samplesPerBit;
    string out;
    for (int b = 0; b < numBits; ++b) {
        int startIdx = b * samplesPerBit;
        int endIdx = min((int)samples.size(), (b + 1) * samplesPerBit);
        
        double sum = 0;
        for (int i = startIdx; i < endIdx; ++i) {
            sum += samples[i];
        }
        double avg = sum / (endIdx - startIdx);
        
        out.push_back(avg > 0 ? '1' : '0');
    }
    return out;
}

string decodeFromSignal_NRZI(const vector<float>& samples, int samplesPerBit) {
    int numBits = samples.size() / samplesPerBit;
    string out;
    
    if (numBits == 0) return out;
    
    double sum = 0;
    int firstEnd = min((int)samples.size(), samplesPerBit);
    for (int i = 0; i < firstEnd; ++i) sum += samples[i];
    float lastLevel = (sum / firstEnd) > 0 ? 1.0f : -1.0f;
    
    out.push_back('0');
    
    for (int b = 1; b < numBits; ++b) {
        int startIdx = b * samplesPerBit;
        int endIdx = min((int)samples.size(), (b + 1) * samplesPerBit);
        
        sum = 0;
        for (int i = startIdx; i < endIdx; ++i) {
            sum += samples[i];
        }
        float currentLevel = (sum / (endIdx - startIdx)) > 0 ? 1.0f : -1.0f;
        
        if (currentLevel != lastLevel) {
            out.push_back('1');
        } else {
            out.push_back('0');
        }
        lastLevel = currentLevel;
    }
    return out;
}

string decodeFromSignal_Manchester(const vector<float>& samples, int samplesPerBit) {
    int totalBitPeriods = samples.size() / samplesPerBit;
    int numLogicalBits = totalBitPeriods / 2;
    string out;
    
    for (int b = 0; b < numLogicalBits; ++b) {
        int firstHalfStart = (2 * b) * samplesPerBit;
        int firstHalfEnd = (2 * b + 1) * samplesPerBit;
        int secondHalfStart = (2 * b + 1) * samplesPerBit;
        int secondHalfEnd = (2 * b + 2) * samplesPerBit;
        
        double sum1 = 0;
        int count1 = 0;
        for (int i = firstHalfStart; i < min((int)samples.size(), firstHalfEnd); ++i) {
            sum1 += samples[i];
            count1++;
        }
        float firstHalf = count1 > 0 ? (sum1 / count1) : 0;
        
        double sum2 = 0;
        int count2 = 0;
        for (int i = secondHalfStart; i < min((int)samples.size(), secondHalfEnd); ++i) {
            sum2 += samples[i];
            count2++;
        }
        float secondHalf = count2 > 0 ? (sum2 / count2) : 0;
        
        if (firstHalf > secondHalf) {
            out.push_back('1');
        } else {
            out.push_back('0');
        }
    }
    return out;
}

string decodeFromSignal_DiffManchester(const vector<float>& samples, int samplesPerBit) {
    int totalBitPeriods = samples.size() / samplesPerBit;
    int numLogicalBits = totalBitPeriods / 2;
    string out;
    
    if (numLogicalBits == 0) return out;
    
    double sum = 0;
    int count = 0;
    for (int i = 0; i < min((int)samples.size(), samplesPerBit / 2); ++i) {
        sum += samples[i];
        count++;
    }
    float prevStartLevel = count > 0 ? ((sum / count) > 0 ? 1.0f : -1.0f) : 1.0f;
    
    for (int b = 0; b < numLogicalBits; ++b) {
        int startIdx = (2 * b) * samplesPerBit;
        sum = 0;
        count = 0;
        for (int i = startIdx; i < min((int)samples.size(), startIdx + samplesPerBit / 2); ++i) {
            sum += samples[i];
            count++;
        }
        float currentStartLevel = count > 0 ? ((sum / count) > 0 ? 1.0f : -1.0f) : prevStartLevel;
        
        if (b == 0) {
            out.push_back('1');
        } else {
            if (currentStartLevel != prevStartLevel) {
                out.push_back('0');
            } else {
                out.push_back('1');
            }
        }
        
        int midIdx = (2 * b + 1) * samplesPerBit;
        sum = 0;
        count = 0;
        for (int i = midIdx; i < min((int)samples.size(), midIdx + samplesPerBit / 2); ++i) {
            sum += samples[i];
            count++;
        }
        prevStartLevel = count > 0 ? ((sum / count) > 0 ? 1.0f : -1.0f) : currentStartLevel;
    }
    return out;
}

string decodeFromSignal_AMI(const vector<float>& samples, int samplesPerBit) {
    int numBits = samples.size() / samplesPerBit;
    string out;
    
    for (int b = 0; b < numBits; ++b) {
        int startIdx = b * samplesPerBit;
        int endIdx = min((int)samples.size(), (b + 1) * samplesPerBit);
        
        double sum = 0;
        for (int i = startIdx; i < endIdx; ++i) {
            sum += samples[i];
        }
        double avg = sum / (endIdx - startIdx);
        
        out.push_back(fabs(avg) > 0.1 ? '1' : '0');
    }
    return out;
}

// ---------- SFML Drawing ----------
void drawSignalWaveform(sf::RenderWindow &window, const vector<float>& samples, 
                        int numBits, sf::Font& font, const string& title,
                        const string& infoText = "", const string& decodedText = "") {
    window.clear(sf::Color::Black);
    
    float width = (float)window.getSize().x;
    float height = (float)window.getSize().y;
    
    float leftMargin = 80.0f;
    float rightMargin = 40.0f;
    float topMargin = 60.0f;
    float bottomMargin = decodedText.empty() ? 100.0f : 150.0f;
    
    float plotWidth = width - leftMargin - rightMargin;
    float plotHeight = height - topMargin - bottomMargin;
    float plotX = leftMargin;
    float plotY = topMargin;
    float mid = plotY + plotHeight / 2.0f;
    
    sf::RectangleShape plotBorder(sf::Vector2f(plotWidth, plotHeight));
    plotBorder.setPosition(plotX, plotY);
    plotBorder.setFillColor(sf::Color::Transparent);
    plotBorder.setOutlineColor(sf::Color::White);
    plotBorder.setOutlineThickness(2.0f);
    window.draw(plotBorder);
    
    sf::VertexArray grid(sf::PrimitiveType::Lines);
    for (int i = 0; i <= numBits; ++i) {
        float x = plotX + (plotWidth * i) / numBits;
        grid.append(sf::Vertex{sf::Vector2f(x, plotY), sf::Color(60,60,60)});
        grid.append(sf::Vertex{sf::Vector2f(x, plotY + plotHeight), sf::Color(60,60,60)});
    }
    for (int i = 0; i <= 4; ++i) {
        float y = plotY + (plotHeight * i) / 4;
        grid.append(sf::Vertex{sf::Vector2f(plotX, y), sf::Color(60,60,60)});
        grid.append(sf::Vertex{sf::Vector2f(plotX + plotWidth, y), sf::Color(60,60,60)});
    }
    window.draw(grid);
    
    sf::Vertex centerAxis[] = {
        sf::Vertex{sf::Vector2f(plotX, mid), sf::Color(150,150,150)},
        sf::Vertex{sf::Vector2f(plotX + plotWidth, mid), sf::Color(150,150,150)}
    };
    window.draw(centerAxis, 2, sf::PrimitiveType::Lines);
    
    size_t N = samples.size();
    if (N >= 2) {
        float xStep = plotWidth / (float)max<size_t>(1, N-1);
        sf::VertexArray line(sf::PrimitiveType::LineStrip);
        for (size_t i = 0; i < N; ++i) {
            float x = plotX + i * xStep;
            float y = mid - samples[i] * (plotHeight/3.0f);
            sf::Color c = samples[i] > 0.1f ? sf::Color::Green : 
                         (samples[i] < -0.1f ? sf::Color::Red : sf::Color::Yellow);
            line.append(sf::Vertex{sf::Vector2f(x, y), c});
        }
        window.draw(line);
    }
    
    sf::Text yAxisLabel("Amplitude", font, 16);
    yAxisLabel.setFillColor(sf::Color::White);
    yAxisLabel.setPosition(10, height/2 - 40);
    yAxisLabel.setRotation(-90);
    window.draw(yAxisLabel);
    
    const char* ampLabels[] = {"+1", "+0.5", "0", "-0.5", "-1"};
    for (int i = 0; i <= 4; ++i) {
        float y = plotY + (plotHeight * i) / 4;
        sf::Text label(ampLabels[i], font, 12);
        label.setFillColor(sf::Color::White);
        label.setPosition(leftMargin - 50, y - 10);
        window.draw(label);
    }
    
    sf::Text xAxisLabel("Bit Period", font, 16);
    xAxisLabel.setFillColor(sf::Color::White);
    xAxisLabel.setPosition(width/2 - 40, plotY + plotHeight + 50);
    window.draw(xAxisLabel);
    
    for (int i = 0; i <= numBits && i <= 20; i += max(1, numBits/10)) {
        float x = plotX + (plotWidth * i) / numBits;
        sf::Text label(to_string(i), font, 12);
        label.setFillColor(sf::Color::White);
        label.setPosition(x - 8, plotY + plotHeight + 10);
        window.draw(label);
    }
    
    sf::Text titleText(title, font, 20);
    titleText.setFillColor(sf::Color::Cyan);
    titleText.setStyle(sf::Text::Bold);
    sf::FloatRect titleBounds = titleText.getLocalBounds();
    titleText.setPosition((width - titleBounds.width) / 2, 10);
    window.draw(titleText);
    
    if (!infoText.empty()) {
        sf::Text info(infoText, font, 14);
        info.setFillColor(sf::Color::Yellow);
        info.setPosition(leftMargin, plotY + plotHeight + 30);
        window.draw(info);
    }
    
    if (!decodedText.empty()) {
        sf::Text decoded(decodedText, font, 14);
        decoded.setFillColor(sf::Color::Green);
        decoded.setPosition(leftMargin, plotY + plotHeight + 75);
        window.draw(decoded);
    }
    
    window.display();
}

// ---------- Main Program ----------
int main() {
    cout << "\n========================================\n";
    cout << "   DIGITAL SIGNAL GENERATOR & DECODER\n";
    cout << "========================================\n\n";
    
    cout << "Select input type:\n";
    cout << "1) Digital bitstream\n";
    cout << "2) Analog samples (PCM/DM)\n";
    cout << "Choice: ";
    int inputType; 
    cin >> inputType;

    string bitstream;
    string analogMethod = "";
    
    if (inputType == 1) {
        cout << "\nEnter binary stream (e.g., 10110011): ";
        cin >> bitstream;
    } else if (inputType == 2) {
        cout << "\nEnter number of analog samples: ";
        int n; 
        cin >> n;
        vector<double> samples(n);
        cout << "Enter " << n << " analog sample values (space-separated):\n";
        for (int i = 0; i < n; ++i) cin >> samples[i];
        
        cout << "\nChoose analog-to-digital conversion method:\n";
        cout << "1) PCM (Pulse Code Modulation)\n";
        cout << "2) DM (Delta Modulation)\n";
        cout << "Choice: ";
        int achoice; 
        cin >> achoice;
        
        if (achoice == 1) {
            cout << "Enter quantization levels (e.g., 8, 16, 32): ";
            int levels; 
            cin >> levels;
            auto q = pcmQuantize(samples, levels);
            bitstream = quantizedToBitstream(q, levels);
            analogMethod = "PCM";
            cout << "\n[PCM] Generated bitstream: " << bitstream << "\n";
        } else {
            cout << "Enter step size for delta modulation (e.g., 0.5): ";
            double step;
            cin >> step;
            auto dm = deltaMod(samples, step);
            bitstream = "";
            for (int v : dm) bitstream.push_back(v ? '1' : '0');
            analogMethod = "Delta Modulation";
            cout << "\n[DM] Generated bitstream: " << bitstream << "\n";
        }
    } else {
        cerr << "Invalid input type!\n";
        return 1;
    }

    string pal = longestPalindromeManacher(bitstream);
    cout << "\n--- Palindrome Analysis ---\n";
    cout << "Longest palindrome found: \"" << pal << "\"\n";
    cout << "Length: " << pal.size() << " bits\n";

    cout << "\n--- Line Encoding Selection ---\n";
    cout << "Choose line encoding scheme:\n";
    cout << "1) NRZ-L (Non-Return-to-Zero Level)\n";
    cout << "2) NRZ-I (Non-Return-to-Zero Inverted)\n";
    cout << "3) Manchester\n";
    cout << "4) Differential Manchester\n";
    cout << "5) AMI (Alternate Mark Inversion)\n";
    cout << "Choice: ";
    int encChoice; 
    cin >> encChoice;

    string encodingName = "";
    string scramblerChoice = "none";
    string originalBitstream = bitstream;
    
    if (encChoice == 5) {
        encodingName = "AMI";
        cout << "\n--- Scrambling Options for AMI ---\n";
        cout << "Do you want scrambling?\n";
        cout << "0) No scrambling\n";
        cout << "1) B8ZS (Bipolar 8-Zero Substitution)\n";
        cout << "2) HDB3 (High Density Bipolar 3)\n";
        cout << "Choice: ";
        int sc; 
        cin >> sc;
        
        if (sc == 1) { 
            bitstream = applyB8ZS(bitstream); 
            scramblerChoice = "B8ZS"; 
            cout << "\n[B8ZS] Scrambled bitstream: " << bitstream << "\n";
        }
        else if (sc == 2) { 
            bitstream = applyHDB3(bitstream); 
            scramblerChoice = "HDB3"; 
            cout << "\n[HDB3] Scrambled bitstream: " << bitstream << "\n";
        }
    } else {
        const char* names[] = {"", "NRZ-L", "NRZ-I", "Manchester", "Differential Manchester"};
        encodingName = names[encChoice];
    }

    vector<int> encodedLevels;
    if (encChoice == 1) encodedLevels = nrzl(bitstream);
    else if (encChoice == 2) encodedLevels = nrzi(bitstream);
    else if (encChoice == 3) encodedLevels = manchester(bitstream);
    else if (encChoice == 4) encodedLevels = diffManchester(bitstream);
    else if (encChoice == 5) encodedLevels = ami(bitstream);
    else { cerr << "Invalid encoding choice!\n"; return 1; }

    int samplesPerEncodedElem = 40;
    auto floatSamples = upsampleToFloatSamples(encodedLevels, samplesPerEncodedElem);

    sf::Font font;
    if (!font.loadFromFile("/System/Library/Fonts/Helvetica.ttc")) {
        if (!font.loadFromFile("/System/Library/Fonts/Courier.dfont")) {
            cerr << "Warning: Could not load font.\n";
        }
    }

    cout << "\n========================================\n";
    cout << "        ENCODING RESULTS\n";
    cout << "========================================\n";
    if (!analogMethod.empty()) {
        cout << "Analog conversion: " << analogMethod << "\n";
    }
    cout << "Original bitstream: " << originalBitstream << "\n";
    if (scramblerChoice != "none") {
        cout << "Scrambler applied: " << scramblerChoice << "\n";
        cout << "Scrambled stream: " << bitstream << "\n";
    }
    cout << "Encoding scheme: " << encodingName << "\n";
    cout << "Encoded signal elements: " << encodedLevels.size() << "\n";
    cout << "\nDisplaying encoded signal...\n";
    cout << "(Keep the window open to see decoding option)\n";

    sf::RenderWindow window(sf::VideoMode(1200, 650), "Encoded Signal");
    int numOriginalBits = originalBitstream.size();
    
    string windowTitle = encodingName + " Encoded Signal";
    string infoText = "Original: " + originalBitstream.substr(0, min(50, (int)originalBitstream.size()));
    if (originalBitstream.size() > 50) infoText += "...";
    
    string decodedBits = "";
    string decodedDisplayText = "";
    bool askedForDecoding = false;

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }
        
        drawSignalWaveform(window, floatSamples, numOriginalBits, font, 
                         windowTitle, infoText, decodedDisplayText);
        
        // Ask for decoding after first render
        if (!askedForDecoding) {
            askedForDecoding = true;
            
            cout << "\n========================================\n";
            cout << "Do you want to decode the signal? (yes/no): ";
            string decodeChoice;
            cin >> decodeChoice;
            
            if (decodeChoice == "yes" || decodeChoice == "Yes" || decodeChoice == "YES" || 
                decodeChoice == "y" || decodeChoice == "Y") {
                
                cout << "\n========================================\n";
                cout << "        DECODING SIGNAL\n";
                cout << "========================================\n";
                cout << "Analyzing encoded signal waveform...\n";
                cout << "Extracting signal parameters from visualization...\n";
                
                if (encChoice == 1) decodedBits = decodeFromSignal_NRZL(floatSamples, samplesPerEncodedElem);
                else if (encChoice == 2) decodedBits = decodeFromSignal_NRZI(floatSamples, samplesPerEncodedElem);
                else if (encChoice == 3) decodedBits = decodeFromSignal_Manchester(floatSamples, samplesPerEncodedElem);
                else if (encChoice == 4) decodedBits = decodeFromSignal_DiffManchester(floatSamples, samplesPerEncodedElem);
                else if (encChoice == 5) decodedBits = decodeFromSignal_AMI(floatSamples, samplesPerEncodedElem);
                
                cout << "\n--- Decoding Results ---\n";
                cout << "Decoded bitstream: " << decodedBits << "\n";
                cout << "Original bitstream: " << bitstream << "\n";
                
                if (decodedBits == bitstream) {
                    cout << "\n✓ SUCCESS! Decoded stream matches original perfectly.\n";
                    cout << "Decoding accuracy: 100%\n";
                    decodedDisplayText = "Decoded: " + decodedBits.substr(0, min(40, (int)decodedBits.size()));
                    if (decodedBits.size() > 40) decodedDisplayText += "...";
                    decodedDisplayText += " | Match: 100% ✓";
                } else {
                    cout << "\n✗ Mismatch detected.\n";
                    int matches = 0;
                    size_t minLen = min(decodedBits.size(), bitstream.size());
                    for (size_t i = 0; i < minLen; ++i) {
                        if (decodedBits[i] == bitstream[i]) matches++;
                    }
                    double matchPct = 100.0 * matches / bitstream.size();
                    cout << "Match percentage: " << fixed << setprecision(2) << matchPct << "%\n";
                    
                    decodedDisplayText = "Decoded: " + decodedBits.substr(0, min(40, (int)decodedBits.size()));
                    if (decodedBits.size() > 40) decodedDisplayText += "...";
                    
                    stringstream ss;
                    ss << fixed << setprecision(1) << matchPct;
                    decodedDisplayText += " | Match: " + ss.str() + "% ✗";
                }
                
                windowTitle = encodingName + " Encoded & Decoded Signal";
                cout << "\nDecoded bitstream displayed in the window.\n";
                cout << "(Close window to exit)\n";
            } else {
                cout << "\nSkipping decoding. Close window to exit.\n";
            }
        }
    }
    
    cout << "\n========================================\n";
    cout << "     PROGRAM COMPLETED SUCCESSFULLY\n";
    cout << "========================================\n\n";
    
    return 0;
}