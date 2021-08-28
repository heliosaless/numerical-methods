#include <iostream>

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

    State operator*(const double s) const{
        return State(v * s,  y * s); 
    }

    friend State operator*(const double s, const State &obj){
        return State(obj.v * s, obj.y * s);
    }
    
    State& operator=(const State &t){
        return *this;
    }

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

    State fim = euler.next(start, delta);
    State med = euler.next(start, delta/2);
    

    return start + delta * ((1./6)*start.derivative(g,m,k) + (4./6)*med.derivative(g,m,k) + (1./6)*fim.derivative(g,m,k));
}

double at(const double time, const double start, const double delta){
    return (time - start)/delta;
}

int main(int argc, char const *argv[])
{
    const double start = 0;
    const double time = 10;
    double delta = 0.1;
    
    RungeKutta rk(10., 2., 1./4);
    State initial(5., 200.);
    State i = initial;
    
    bool cond1 = false;
    bool cond2 = false;
    int stateNumber = 0;

    while (delta >= 0.0001)
    {
        std::cout << "Usando delta = " << delta << std::endl;

        int count = at(time, start, delta);
        State i2 = initial;

        while (--count >= 0)
        {
            State a = rk.next(i2, delta);
            if (count == 0)
            {
                std::cout << "\tNo tempo 10 temos valor " << a.getY() << std::endl;
            }
            i2.setV(a.getV());
            i2.setY(a.getY());
        }

        delta = delta / 10;
    }

    std::cout << std::endl;
    delta = 0.1;
    while(true){
        ++stateNumber;
        State a = rk.next(i, delta);

        if(!cond1 && a.getY() < i.getY()){
            std::cout << "\tA altura maxima alcancada foi " << i.getY() << "m";
            std::cout << " Essa altura ocorreu " << (stateNumber-1)*delta << "s apos o lancamento." << std::endl;
            cond1 = true;
        }

        if(!cond2 && a.getY() < 0){
            double tempo = ((start + (stateNumber-1)*delta)+(start + stateNumber*delta))/2;
            std::cout << "\tCruzou o mar no tempo " << tempo << "s";
            std::cout << " com velocidade de " << (i.getV() + a.getV())/2 << "m/s." << std::endl;
            cond2 = true;
        }

        if(cond1 && cond2) break;
        i.setV(a.getV());
        i.setY(a.getY());
    }

    

    return 0;
}
