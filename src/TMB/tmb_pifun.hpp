
//Removal pi function with constant time intervals
template<class Type>
vector<Type> pifun_removal(vector<Type> p){
  int J = p.size();
  vector<Type> pi(J);
  pi(0) = p(0);
  for (int j=1; j<J; j++){
    pi(j) = pi(j-1) / p(j-1) * (1-p(j-1)) * p(j);
  }
  return(pi);
}

//Removal pi function with variable time intervals
template<class Type>
vector<Type> pifun_removal(vector<Type> p, vector<int> times){
  int J = p.size();
  for (int j=0; j<J; j++){
    p(j) = 1 - pow(1 - p(j), times(j));
  }
  return pifun_removal(p);
}

//Double observer pi function
//p *must* be length 2
template<class Type>
vector<Type> pifun_double(vector<Type> p){
  vector<Type> pi(3);
  pi(0) = p(0) * (1 - p(1));
  pi(1) = p(1) * (1 - p(0));
  pi(2) = p(0) * p(1);
  return(pi);
}

//Dependent double observer pi function
//p *must* be length 2
template<class Type>
vector<Type> pifun_dep_double(vector<Type> p){
  vector<Type> pi(2);
  pi(0) = p(0);
  pi(1) = p(1) * (1 - p(0));
  return(pi);
}

//Umbrella function
template<class Type>
vector<Type> pifun(vector<Type> p, int pifun_type){
  if(pifun_type == 0){
    return(pifun_removal(p));
  } else if(pifun_type == 1){
    return(pifun_double(p));
  } else if(pifun_type == 2){
    return(pifun_dep_double(p));
  } else {
    throw(std::invalid_argument("invalid pifun"));
  }
}
