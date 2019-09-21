//
// Created by Lander on 5/8/18.
//

#include "session1.cpp"

std::string replaceChar3D(const LParser::LSystem3D &l_system, const char c, unsigned int amount) {
    if (c == '+' || c == '-' || c == '^' || c == '&' || c == '\\' || c == '/' || c == '|' || c == '(' || c == ')') {
        std::string w;
        return w + c;
    }
    if (amount == 1) {
        return l_system.get_replacement(c);
    }
    std::string replaced;
    for (char w : l_system.get_replacement(c)) {
        replaced.append(replaceChar3D(l_system, w, amount - 1));
    }
    return replaced;
}

void generateLSystem3D(Figure &ImageFigure, const LParser::LSystem3D &l_system) {
    std::string drawstring;
    Face g;
    std::stack<Vector3D> pos;
    std::stack<Vector3D> sH;
    std::stack<Vector3D> sL;
    std::stack<Vector3D> sU;
    double angle = convertToRadiant(l_system.get_angle());
    Vector3D H;
    H.x = 1;
    H.y = 0;
    H.z = 0;
    Vector3D Hn;
    Vector3D L;
    L.x = 0;
    L.y = 1;
    L.z = 0;
    Vector3D Ln;
    Vector3D U;
    U.x = 0;
    U.y = 0;
    U.z = 1;
    Vector3D Un;
    Vector3D current;
    current.x = 0;
    current.y = 0;
    current.z = 0;
    Vector3D next;

    if (l_system.get_nr_iterations() == 0) {
        drawstring = l_system.get_initiator();
    } else {
        for (char c : l_system.get_initiator()) {
            drawstring.append(replaceChar3D(l_system, c, l_system.get_nr_iterations()));
        }
    }

    for (char c : drawstring) {
        if (c == '+') {
            Hn = H * std::cos(angle) + L * std::sin(angle);
            Ln = -H * std::sin(angle) + L * std::cos(angle);
            H = Hn;
            L = Ln;
            continue;
        }
        if (c == '-') {
            Hn = H * std::cos(-angle) + L * std::sin(-angle);
            Ln = -H * std::sin(-angle) + L * std::cos(-angle);
            H = Hn;
            L = Ln;
            continue;
        }
        if (c == '^') {
            Hn = H * std::cos(angle) + U * std::sin(angle);
            Un = -H * std::sin(angle) + U * std::cos(angle);
            H = Hn;
            U = Un;
            continue;
        }
        if (c == '&') {
            Hn = H * std::cos(-angle) + U * std::sin(-angle);
            Un = -H * std::sin(-angle) + U * std::cos(-angle);
            H = Hn;
            U = Un;
            continue;
        }
        if (c == '\\') {
            Ln = L * std::cos(angle) - U * std::sin(angle);
            Un = L * std::sin(angle) + U * std::cos(angle);
            L = Ln;
            U = Un;
            continue;
        }
        if (c == '/') {
            Ln = L * std::cos(-angle) - U * std::sin(-angle);
            Un = L * std::sin(-angle) + U * std::cos(-angle);
            L = Ln;
            U = Un;
            continue;
        }
        if (c == '|') {
            H = -H;
            L = -L;
            continue;
        }
        if (c == '(') {
            pos.push(current);
            sH.push(H);
            sU.push(U);
            sL.push(L);
            continue;
        }
        if (c == ')') {
            current = pos.top();
            pos.pop();
            H = sH.top();
            sH.pop();
            U = sU.top();
            sU.pop();
            L = sL.top();
            sL.pop();
            continue;
        }

        next = current + H;

        if (l_system.draw(c)) {
            ImageFigure.points.push_back(current);
            ImageFigure.points.push_back(next);
            g.index = {ImageFigure.points.size() - 2, ImageFigure.points.size() - 1};
            ImageFigure.faces.push_back(g);
        }

        current = next;
    }
};
