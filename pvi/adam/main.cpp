#include <iostream>
#include <math.h>
class State
{
private:
    double v, y;

public:
    State(const double v, const double y);
    State();
    State derivative(const double g, const double m, const double k);

    State operator+(State const &obj){
        State res;
        res.v = v + obj.v;
        res.y = y + obj.y;
        return res;
    }

    State operator-(State const &obj){
        State res;
        res.v = v - obj.v;
        res.y = y - obj.y;
        return res;
    }

    State operator/(const State &obj) const{
        State res;
        res.v = v / obj.v;
        res.y = y / obj.y;
        return res;
    }

    State operator*(const double s) const{
        return State(v * s,  y * s); 
    }

    friend State operator*(const double s, const State &obj){
        return State(obj.v * s, obj.y * s);
    }
    
    State& operator=(const State &t){
        return *this;
    }

    const double mod() {return sqrt(v*v+y*y);}
    const double getV(){ return v; }
    const double getY(){ return y; }
    void setV(const double v_) { v = v_; }
    void setY(const double y_) { y = y_; }
};

State::State(){v = 0; y = 0;}
State::State(const double v, const double y):v(v),y(y){}
State State::derivative(const double g, const double m, const double k) { return State(-g - (k / m) * v, v); }

class Euler
{
private:
    double g, m, k;
public:
    Euler() {}
    Euler(const double g, const double m, const double k);
    State next(State start, const double delta);
};

Euler::Euler(const double g, const double m, const double k):g(g),m(m),k(k){}
State Euler::next(State start, const double delta){
    return start + delta*start.derivative(g, m, k);
}

class RungeKutta
{
private:
    double g, m, k;
    Euler euler;
public:
    RungeKutta(const double g, const double m, const double k);
    State next(State start, const double delta);
};

RungeKutta::RungeKutta(const double g, const double m, const double k) : g(g), m(m), k(k) {euler = Euler (g,m,k);}
State RungeKutta::next(State start, const double delta){

    State s2 = start + (delta/2)*start.derivative(g,m,k);
    State s3 = start + (delta/2)*s2.derivative(g,m,k);
    State s4 = start + (delta)*s3.derivative(g,m,k);
    
    return start + (delta/6) * (start.derivative(g,m,k) + 2*s2.derivative(g,m,k) + 2*s3.derivative(g,m,k) + s4.derivative(g,m,k));
}

double at(const double time, const double start, const double delta){
    return (time - start)/delta;
}

class AdamsMoulton{
private:
    const double g,m,k;
public:
    AdamsMoulton(const double g, const double m, const double k);
    State next(State s3, const double delta, State s2, State s1, State s0, const double epls);
};

AdamsMoulton::AdamsMoulton(const double g, const double m, const double k) : g(g), m(m), k(k) {}
State AdamsMoulton::next(State s3, const double delta, State s2, State s1, State s0, const double epls = 10e-6){

    State aux = s3 + (delta/24) * ((55.)*s3.derivative(g,m,k) + (-59.)*
                                s2.derivative(g,m,k) + (37.)*s1.derivative(g,m,k) + (-9.)*s0.derivative(g,m,k));

    double err = 1;
    while(true){
        State a = s3 + (delta/24)*((9.)*aux.derivative(g,m,k) + 19*s3.derivative(g,m,k)
                                             + (-5)*s2.derivative(g,m,k) + s1.derivative(g,m,k));

        
        if(((a - aux)/a).mod() < epls) return a;

        aux.setV(a.getV());
        aux.setY(a.getY());
    }

    return State();
}



int main(int argc, char const *argv[])
{
    const double start = 0;
    const double time = 10;
    double delta = 0.1;
    
    RungeKutta rk(10., 2., 1./4);
    AdamsMoulton adam(10., 2., 1./4);
    State initial(5., 200.);

    while (delta >= 0.0001)
    {
        std::cout << "Usando delta = " << delta << std::endl;

        int count = at(time, start, delta);
        State s0 = initial;
        State s1 = rk.next(s0, delta);
        State s2 = rk.next(s1, delta);
        State s3 = rk.next(s2, delta);

        while (--count >= 0)
        {
            State a = adam.next(s3, delta, s2, s1, s0);
            if (count == 0)
                std::cout << "\tNo tempo 10 temos valor " << a.getY() << std::endl;
            
            s0.setV(s1.getV());
            s0.setY(s1.getY());
            s1.setV(s2.getV());
            s1.setY(s2.getY());
            s2.setV(s3.getV());
            s2.setY(s3.getY());
            s3.setV(a.getV());
            s3.setY(a.getY());
        }

        delta = delta / 10;
    }

    delta = 0.1;
    int stateNumber = 3;
    bool cond1 = false;
    bool cond2 = false;

    State s0 = initial;
    State s1 = rk.next(s0, delta);
    State s2 = rk.next(s1, delta);
    State s3 = rk.next(s2, delta);

    std::cout << std::endl;

    while(true){
        ++stateNumber;
        State a = adam.next(s3, delta, s2, s1, s0);

        if(!cond1 && a.getY() < s3.getY()){
            std::cout << "\tA altura maxima alcancada foi " << s3.getY() << "m";
            std::cout << " Essa altura ocorreu " << (stateNumber-1)*delta << "s apos o lancamento." << std::endl;
            cond1 = true;
        }

        if(!cond2 && a.getY() < 0){
            double tempo = ((start + (stateNumber-1)*delta)+(start + stateNumber*delta))/2;
            std::cout << "\tCruzou o mar no tempo " << tempo << "s";
            std::cout << " com velocidade de " << (s3.getV() + a.getV())/2 << "m/s." << std::endl;
            cond2 = true;
        }

        if(cond1 && cond2) break;
        s0.setV(s1.getV());
        s0.setY(s1.getY());
        s1.setV(s2.getV());
        s1.setY(s2.getY());
        s2.setV(s3.getV());
        s2.setY(s3.getY());
        s3.setV(a.getV());
        s3.setY(a.getY());
    }

    return 0;
}
