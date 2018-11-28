pipeline {
  agent {
    docker {
      image 'python:3.7-slim'
    }

  }
  stages {
    stage('') {
      steps {
        sh 'python -c "import uuid;print(uuid.uuid4());"'
      }
    }
  }
}