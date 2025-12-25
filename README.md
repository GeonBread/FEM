# FEM
Finite Element Method Code 

문제 및 과제 풀이할 때 작성한 MATLAB Code인데 그냥 버리기 아까워서 올려둚

book : A First Course in the Finite Element Method

2025-12-03
Grid Element에 대한 문제 풀이 코드 업로드 함. 250개에 251개의 노드를 가진 복잡한 Grid Element에 대한 해석. 교수님이 작성하신 코드를 바탕으로 풀이.
기본적으로 교수님 풀이는 Element의 Node Number와 좌표 부여 시 불규칙성에 의해 풀이 실수를 할 가능성이 큼. 하지만 내가 만든 방식대로 하면 250개 정도의 element도
쉽게 코드로 부여할 수 있게 됨. 다만 좌표는 일일이 해야 함.

2025-12-26
종강하고 업로드.

Frame Element code인데 Element 나눠서 쉽게 matE, matP 짜는 방식으로 풀었음. 
Chapter 6부터는 2D Element에 관하여 다루고, 학부 수준에서는 Linear System만 다루기 때문에 Constant Strain Triangle인 CST Element와
Ractangle Element인 Q4 Element에 대해서 다뤘음.

CST, Q4 feature?
1. displacement에 대한 interpolation이 linear 식임. 따라서 당연히 non-linear한 부분을 다룰 수 없고, 실제 거동과 너무 다르기 때문에 한계가 엄청 큼.
2. 하지만 선형 element라서 배우기는 쉬운 편임.
3. 중요한 건 여기서부턴 Equilibrium으로 Local Stiffness Equation을 derive 할 수 없기 때문에 Minimum Potential Energy Theory를 사용함.

2D Element 다루면서 확실히 퍼텐셜 에너지라는 것에 더 익숙해진 느낌임. 에너지 관점에서 뭔가를 해석한다는 건 정말 쉽고 간단한 작업이 될 수 있다는 것도 느꼈음.

FEM에 대하여?
- 단순히 수치 풀이 방식인 유한요소해석을 배우는 것을 넘어서 고체역학적 지식을 다시 한번 복기할 수 있는 기회가 된 것 같음.
- 구체적으로 학부 과정 때 배우던 모델들(classical Torsion model, Euler beam)에 대해서도 다시 한번 상기할 수 있었고, 유도 과정을 보면서 평형식 세우는 방법에 대해 더 잘 알게 된 것 같음.
항상 헷갈리던게 평형식 부호 문제였는데, 이 점에 대해서도 좀 더 잘 알게 된 것 같음. 또 MATLAB coding이 이젠 완전 익숙해짐. 디버깅도 visual studio code로 하던 거랑 다를 게 없음.
무엇보다 python의 numpy가 MATLAB 라이브러리를 완전히 배낀 거였다니.. 정말 똑같다.

강의 후기?
- 기계과 수업 중 가장 어려운 수업 중 하나로 뽑히는 수업인데, 맞긴 한 듯. 다만 난 코딩을 원래부터 할 줄 알기도 했고, 솔직히 이번 학기 전공 수업들 중 이게 가장 재밌어서 할만했음.
중간고사 때는 근데 벽 느끼긴 함. 너무 어렵고 시험 시간이 4시간이나 되니까 힘들었음. 그래도 FEM은 전산 유체 역학과 함께 꼭 듣고 싶었던 과목 중 하나이기도 했음. 왜냐하면 맨날
2학년 때 손으로 풀어서 Displacement 구하고 Stress구하던 걸 실제 프로그램 돌려보면서 구체적으로 어떻게 변하는지도 보고 싶었기 때문임. 결과적으로 아주 만족스러운 수업 중 하나.
또 저학년 때 대충 넘어간 개념들에 대해서도 다시 한번 상기하는 걸 넘어서 더 깊게 이해할 수 있는 시간을 가진 것 같음. 역시 이래서 공부는 끝이 없고, 같은 책을 여러 번 봐도
항상 배우는 게 있다더니 맞는 말인 것 같음. 그동안 안다고 생각했던 건 사실 그냥 외웠던 거였음. 공부에 대한 자세도 좀 배운 듯. 
